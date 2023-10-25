function c_alphas_to_points(chain::BiochemicalAlgorithms.Chain)
    c_alphas = []
    point_to_residue_indices = []
    for (i, fragment) in enumerate(fragments(chain))
        if(is_amino_acid(fragment.name))
            ca = filter(x -> x.element==Elements.C && x.name=="CA", atoms(fragment))
            @assert length(ca)==1
            push!(c_alphas, ca[1].r)
            push!(point_to_residue_indices, i)
        end
    end
    sphere_radius = 0.05
    sphere_mesh = simplexify(Sphere{3, Float32}((0,0,0), sphere_radius))
    spheres = map(a -> Translate(Float32.(a)...)((sphere_mesh)), c_alphas)
    spheres = map(s -> ColoredMesh(s, (0, 0, 0)), spheres)
    m1 = reduce(merge, spheres)
    export_mesh_to_ply("c-alpha.ply", m1)

    return c_alphas, point_to_residue_indices
end

function generate_points_carson_bugg(chain::BiochemicalAlgorithms.Chain, offset_helix_points::Bool)
    # find relevant atom positions
    c_alphas = []
    oxygens = []
    structures = []
    point_to_residue_indices = []
    for (i, fragment) in enumerate(fragments(chain))
        if(is_amino_acid(fragment.name))
            ca = filter(x -> x.element==Elements.C && x.name=="CA", atoms(fragment))
            ox = filter(x -> x.element==Elements.O && x.name=="O", atoms(fragment))
            @assert length(ca)==1 && length(ox)==1

            push!(c_alphas, ca[1])
            push!(oxygens, ox[1])
            push!(structures, fragment.properties[:SS])
            push!(point_to_residue_indices, i)
        end
    end

    sphere_radius = 0.05
    sphere_mesh = simplexify(Sphere{3, Float32}((0,0,0), sphere_radius))
    spheres = map(a -> Translate(Float32.(a.r)...)((sphere_mesh)), c_alphas)
    spheres = map(s -> ColoredMesh(s, (0, 0, 0)), spheres)
    m1 = reduce(merge, spheres)
    export_mesh_to_ply("c-alpha.ply", m1)

    main_points = Matrix{Float64}(undef, 3, length(c_alphas)-1)
    minor_points = Matrix{Float64}(undef, 3, length(c_alphas)-1)

    # construct planes to obtain spline control points
    current_flip = false
    prev_D = nothing
    for i=1:length(c_alphas)-1
        A = c_alphas[i+1].r - c_alphas[i].r
        B = oxygens[i].r - c_alphas[i].r
        C = cross(A, B)
        D = cross(C, A)

        C = C/norm(C)
        D = D/norm(D)

        flip = prev_D!==nothing && abs(acos(dot(prev_D, D)))>0.5*π
        if(flip)
            current_flip = !current_flip
        end
        prev_D = D
        
        P = c_alphas[i].r + 0.5 * A
        if(offset_helix_points && (structures[i]==BiochemicalAlgorithms.SecondaryStructure.HELIX || structures[i+1]==BiochemicalAlgorithms.SecondaryStructure.HELIX))
            P += 1.5 * C
        end

        main_points[:, i] = P
        if(current_flip)
            minor_points[:, i] = P+D
        else
            minor_points[:, i] = P-D
        end
    end
    sphere_radius = 0.1
    sphere_mesh = simplexify(Sphere{3, Float32}((0,0,0), sphere_radius))
    spheres = map(a -> Translate(Float32.(a)...)((sphere_mesh)), eachcol(main_points))
    spheres = map(s -> ColoredMesh(s, (0, 0, 255)), spheres)
    m1 = reduce(merge, spheres)
    export_mesh_to_ply("main_points.ply", m1)
    return main_points, minor_points, point_to_residue_indices
end

function num_points(spline, resolution)
    dict_key = Int(round(resolution*1000))
    if(dict_key ∈ keys(spline.num_points_per_resolution))
        return spline.num_points_per_resolution[dict_key]
    else
        num_points = []
        sample_mapping = []
        i = 1
        while i+3 <= size(spline.controlPoints, 2)
            distance = norm(spline.controlPoints[:, i+1] .- spline.controlPoints[:, i+2])
            push!(num_points, max(2, convert(Int, ceil(resolution * 1 * distance)))) #TODO factor

            if(spline.controlPointStrategy==ControlPoints.C_ALPHA)
                first_half_num = num_points[end] ÷ 2
                second_half_num = num_points[end]-first_half_num
                if(i+3!=size(spline.controlPoints, 2))
                    second_half_num -= 1
                end

                push!(sample_mapping, 
                        repeat([spline.point_to_residue_indices[i+1]], first_half_num)..., 
                        repeat([spline.point_to_residue_indices[i+2]], second_half_num)...)
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
        spline.num_points_per_resolution[dict_key] = num_points
        spline.sample_mapping_per_resolution[dict_key] = sample_mapping
        return num_points
    end
end

function evaluate_generic_quadruple_spline(control_points, num_points, fn)
    result_points = Matrix(undef, 3, sum(num_points)-length(num_points)+1)
    i = 1
    a = 1
    while i+3 <= size(control_points, 2) # loop over quadruples of control_points
        points = fn((control_points[:, i], control_points[:, i+1], control_points[:, i+2], control_points[:, i+3]), num_points[i])
        result_points[:, a:a+num_points[i]-1] = points
        a += num_points[i]-1
        i += 1  
    end

    return result_points
end