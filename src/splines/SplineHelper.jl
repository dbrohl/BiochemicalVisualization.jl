function get_c_alpha_positions(chain::Chain{T}, with_oxygens=false) where T
    cas::Vector{Tuple{Atom{T}, Fragment{T}}} = []
    os ::Vector{Tuple{Atom{T}, Fragment{T}}}= []
    for x in eachatom(chain)
        if x.element==Elements.C && x.name=="CA"
            parent = parent_fragment(x)
            if is_amino_acid(parent)
                push!(cas, (x, parent))
            end
        end
    
        if with_oxygens && x.element==Elements.O && x.name=="O"
            parent = parent_fragment(x)
            if is_amino_acid(parent)
                push!(os, (x, parent))
            end
        end
    end

    if(with_oxygens && length(cas)!=length(os))
        log_warning("different number of c-alpha ($(length(cas))) and oxygen ($(length(os))) atoms")
        if(length(cas)<length(os))
            os = os[1:length(cas)]
        else
            cas = cas[1:length(os)]
        end
    end

    positions = Matrix{T}(undef, 3, length(cas))
    indices = Vector{Int}(undef, length(cas))
    residue_info_dict = Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}}()

    for (i, (ca, fragment)) in enumerate(cas)
        positions[:, i] = ca.r
        indices[i] = fragment.idx
        residue_info_dict[fragment.idx] = (fragment.name, fragment.properties[:SS])
    end
    if(with_oxygens)
        oxygen_positions = Matrix{T}(undef, 3, length(os))
        for (i, (o, fragment)) in enumerate(os)
            oxygen_positions[:, i] = o.r
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

        C ./= norm(C)
        D ./= norm(D)

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

function calculate_resolution_dependent_data(spline, resolution)
    num_points::Vector{Int} = []
    sample_mapping::Vector{Int} = []
    i = 1
    while i+3 <= size(spline.controlPoints, 2)
        distance = @views norm(spline.controlPoints[:, i+1] .- spline.controlPoints[:, i+2])
        push!(num_points, max(2, convert(Int, ceil(resolution * 1 * distance)))) #TODO factor

        if(spline.controlPointStrategy==ControlPoints.C_ALPHA)
            first_half_num = num_points[end] ÷ 2
            second_half_num = num_points[end]-first_half_num
            if(i+3!=size(spline.controlPoints, 2))
                second_half_num -= 1
            end

            push!(sample_mapping, 
                    repeat([spline.point_to_residue_indices[i+1]], first_half_num)..., 
                    repeat([spline.point_to_residue_indices[i+2]], second_half_num)...) #alloc
        elseif(spline.controlPointStrategy==ControlPoints.MID_POINTS)
            repeats = num_points[end]
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