function get_c_alpha_positions(chain::Chain{T}, with_oxygens=false) where T

    sys = parent_system(chain)

    atom_data = atoms_df(chain)[:, [:idx, :element, :name, :r, :properties]]
    subset!(atom_data, :name => ByRow(n -> n=="CA" || (with_oxygens && n=="O")))
    select!(atom_data, :idx, :element, :name, :r, :properties, :idx => ByRow(id -> parent_fragment(atom_by_idx(sys, id)).idx) => :parent_idx)

    residue_data = fragments_df(chain)[:, [:idx, :name, :properties]]
    subset!(residue_data, :name => names -> is_amino_acid.(names))
    
    result_data = innerjoin(atom_data, residue_data, on = :parent_idx => :idx, renamecols = "_atom" => "_residue")

    if with_oxygens
        n = size(result_data, 1) ÷ 2
        oxygen_positions = Matrix{T}(undef, 3, n)
    else
        n = size(result_data, 1)
    end
    positions = Matrix{T}(undef, 3, n)
    indices = Vector{Int}(undef, n)
    residue_info_dict = Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}}()

    c_count = 1
    o_count = 1
    for row in eachrow(result_data)
        if row.element_atom==Elements.C
            positions[:, c_count] = row.r_atom
            indices[c_count] = row.parent_idx
            residue_info_dict[row.parent_idx] = (row.name_residue, row.properties_residue[:SS])
            c_count += 1
        elseif with_oxygens
            oxygen_positions[:, o_count] = row.r_atom
            o_count += 1
        end
    end

    if(!issorted(indices))
        order = sortperm(indices)
        indices = indices[order]
        positions = positions[:, order]
        if(with_oxygens)
            oxygen_positions = oxygen_positions[:, order]
        end
    end

    if(with_oxygens)
        return positions, indices, residue_info_dict, oxygen_positions
    else 
        return positions, indices, residue_info_dict
    end
end

function generate_points_carson_bugg(chain::Chain{T}, offset_helix_points::Bool) where T
    # find relevant atom positions
    c_positions, point_to_residue_indices, residue_info_dict, o_positions = get_c_alpha_positions(chain, true)
    
    structures = Vector{BiochemicalAlgorithms.SecondaryStructure.T}(undef, length(point_to_residue_indices))
    for (i, res_idx) in enumerate(point_to_residue_indices)
        structures[i] = residue_info_dict[res_idx][2]
    end
    if(size(c_positions, 2)<3)
        throw(ErrorException("too few ($(size(c_positions, 2))) c_alpha atoms to compute spline with Carson&Bugg method"))
    end

    # sphere_radius = 0.05
    # sphere_mesh = simplexify(Sphere{3, Float32}((0,0,0), sphere_radius))
    # spheres = map(a -> Translate(Float32.(a.r)...)((sphere_mesh)), c_alphas)
    # spheres = map(s -> ColoredMesh(s, (0, 0, 0)), spheres)
    # m1 = reduce(merge, spheres)
    # export_mesh_to_ply("c-alpha.ply", m1)

    main_points = Matrix{T}(undef, 3, size(c_positions, 2)-1)
    minor_points = Matrix{T}(undef, 3, size(c_positions, 2)-1)

    # construct planes to obtain spline control points
    current_flip = false
    prev_D::Union{Vector{T}, Nothing} = nothing
    for i=1:size(c_positions, 2)-1
        A = c_positions[:, i+1] .- c_positions[:, i]
        B = o_positions[:, i] .- c_positions[:, i]
        C = cross(A, B)
        D = cross(C, A)

        normalize!(C)
        normalize!(D)

        flip = prev_D!==nothing && abs(acos(dot(prev_D, D)))>0.5*π
        if(flip)
            current_flip = !current_flip
        end
        prev_D = D
        
        P = c_positions[:,i] + T(0.5) * A
        if(offset_helix_points && (structures[i]==BiochemicalAlgorithms.SecondaryStructure.HELIX || structures[i+1]==BiochemicalAlgorithms.SecondaryStructure.HELIX))
            P += T(1.5) * C
        end

        main_points[:, i] = P
        if(current_flip)
            minor_points[:, i] = P+D
        else
            minor_points[:, i] = P-D
        end
    end
    # sphere_radius = 0.1
    # sphere_mesh = simplexify(Sphere{3, Float32}((0,0,0), sphere_radius))
    # spheres = map(a -> Translate(Float32.(a)...)((sphere_mesh)), eachcol(main_points))
    # spheres = map(s -> ColoredMesh(s, (0, 0, 255)), spheres)
    # m1 = reduce(merge, spheres)
    # export_mesh_to_ply("main_points.ply", m1)
    return main_points, minor_points, point_to_residue_indices, residue_info_dict
end

function calculate_resolution_dependent_data(spline::Union{CatmullRom, CubicB}, resolution)
    num_points = Vector{Int}(undef, size(spline.controlPoints, 2)-3)
    sample_mapping::Vector{Int} = []

    sizehint!(sample_mapping, Int(round((size(spline.controlPoints, 2)-3) * 3.75 * resolution)))
    i = 1
    while i+3 <= size(spline.controlPoints, 2)
        distance = @views norm(spline.controlPoints[:, i+1] .- spline.controlPoints[:, i+2])
        num_points[i] = max(2, convert(Int, ceil(resolution * distance)))

        if(spline.controlPointStrategy==ControlPoints.C_ALPHA)
            first_half_num = num_points[i] ÷ 2
            second_half_num = num_points[i]-first_half_num
            if(i+3!=size(spline.controlPoints, 2))
                second_half_num -= 1
            end

            push!(sample_mapping, 
                    repeat([spline.point_to_residue_indices[i+1]], first_half_num)..., 
                    repeat([spline.point_to_residue_indices[i+2]], second_half_num)...) #alloc
        elseif(spline.controlPointStrategy==ControlPoints.MID_POINTS)
            repeats = num_points[i]
            if(i+3!=size(spline.controlPoints, 2))
                repeats -= 1
            end

            push!(sample_mapping, 
                    repeat([spline.point_to_residue_indices[i+2]], repeats)...)
            
        end
        i += 1
        
    end

    dict_key = Int(round(resolution*1000))
    spline.num_points_per_resolution[dict_key] = num_points
    spline.sample_mapping_per_resolution[dict_key] = sample_mapping
end

function num_points(spline, resolution)
    dict_key = Int(round(resolution*1000))
    if(dict_key ∉ keys(spline.num_points_per_resolution))
        calculate_resolution_dependent_data(spline, resolution)
    end
    return spline.num_points_per_resolution[dict_key]
end

function sample_to_fragment_index_mapping(spline, resolution)
    dict_key = Int(round(resolution*1000))
    if(dict_key ∉ keys(spline.sample_mapping_per_resolution))
        calculate_resolution_dependent_data(spline, resolution)
    end
    return spline.sample_mapping_per_resolution[dict_key]
end

function evaluate_generic_quadruple_spline(control_points::Matrix{T}, num_points, fn) where T
    result_points = Matrix{T}(undef, 3, sum(num_points)-length(num_points)+1)
    i = 1
    a = 1
    while i+3 <= size(control_points, 2) # loop over quadruples of control_points
        result_points[:, a:a+num_points[i]-1] = @views fn((control_points[:, i], control_points[:, i+1], control_points[:, i+2], control_points[:, i+3]), num_points[i])
        a += num_points[i]-1
        i += 1  
    end

    return result_points
end