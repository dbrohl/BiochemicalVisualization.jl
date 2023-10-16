mutable struct CubicB
    main_ribbon::AbstractMatrix
    minor_ribbon_a::AbstractMatrix
    minor_ribbon_b::AbstractMatrix

    function CubicB(chain)
        ribbon_width = Dict(SecondaryStructure.HELIX => 3, SecondaryStructure.SHEET => 3, SecondaryStructure.NONE => 1)
    
        c_alphas = []
        oxygens = []
        for fragment in eachfragment(chain)
            if(is_amino_acid(fragment.name))
                ca = filter(x -> x.element==Elements.C && x.name=="CA", atoms(fragment))
                ox = filter(x -> x.element==Elements.O && x.name=="O", atoms(fragment))
                @assert length(ca)==1 && length(ox)==1

                push!(c_alphas, ca[1])
                push!(oxygens, ox[1])
            end
        end
        structures = [SecondaryStructure.NONE for a in c_alphas] # TODO 

        # template for debug output
        # colors = [(255, 0, 0), (255, 255, 0), (0, 255, 0), (0, 255, 255), (0, 0, 255), (255, 0, 255)]
        sphere_radius = 0.2
        sphere_mesh = discretize(Sphere{3, Float64}((0,0,0), sphere_radius), RegularDiscretization(6))
        # debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a.r)...)((sphere_mesh)), (0, 0, 0)), c_alphas))
        # export_mesh_to_ply("c_alphas.ply", debug_mesh)
        # debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a.r)...)((sphere_mesh)), (0, 0, 255)), oxygens))
        # export_mesh_to_ply("oxygens.ply", debug_mesh)

        main_points = Matrix{Float64}(undef, 3, length(c_alphas)-1)
        outer_points = Matrix{Float64}(undef, 3, length(c_alphas)-1)
        outer_points_2 = Matrix{Float64}(undef, 3, length(c_alphas)-1)
        prev_D = nothing

        current_flip = false
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
            # if(prev_D !== nothing)
            #     angle = acos(dot(prev_D, D))
            #     println(angle/(2*π)*360, " ", abs(angle)<=0.5*π)
            # end
            prev_D = D
            
            P = c_alphas[i].r + 0.5 * A
            if(structures[i]==SecondaryStructure.HELIX || structures[i+1]==SecondaryStructure.HELIX) # TODO correct condition?
                # TODO translate
            end
            D *= 0.5 * ribbon_width[structures[i]]
    
            main_points[:, i] = P
            if(current_flip)
                outer_points[:, i] = P+D
                outer_points_2[:, i] = P-D
            else
                outer_points[:, i] = P-D
                outer_points_2[:, i] = P+D
            end 

            
        end

        dm1 = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (100, 100, 100)), eachcol(main_points)))
        dm2 = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (255, 0, 0)), eachcol(outer_points)))
        dm3 = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (0, 255, 0)), eachcol(outer_points_2)))
        debug_mesh = BiochemicalVisualization.merge(BiochemicalVisualization.merge(dm1, dm2), dm3)
        export_mesh_to_ply("control_points.ply", debug_mesh)

        new(main_points,outer_points, outer_points_2)
    end
end

function (spline::CubicB)(resolution; with_frames=false)

    num_points = []
    i = 1
    while i+3 <= size(spline.main_ribbon, 2)
        distance = norm(spline.main_ribbon[:, i+1] .- spline.main_ribbon[:, i+2])
        push!(num_points, max(2, convert(Int, ceil(resolution * 2 * distance))))
        i += 1
    end

    points = Matrix(undef, 3, sum(num_points)-length(num_points)+1)
    velocities = Matrix(undef, 3, sum(num_points)-length(num_points)+1)
    accelerations = Matrix(undef, 3, sum(num_points)-length(num_points)+1)
    if(with_frames)
        rs = Matrix(undef, 3, sum(num_points)-length(num_points)+1) # points of minor ribbon
    end

    i = 1
    a = 1
    M = [1 4 1 0
        -3 0 3 0
        3 -6 3 0
        -1 3 -3 1]
    while i+3 <= size(spline.main_ribbon, 2) # loop over quadruples of controlPoints

        p_matrix = spline.main_ribbon[:, i:i+3]'
        sampling_range = range(0, 1, num_points[i])
        for t in sampling_range
            points[:, a] = [1 t t^2 t^3] * 1/6 * M * p_matrix
            velocities[:, a] = [0 1 2*t 3*t^2] * 1/6 * M * p_matrix
            accelerations[:, a] = [0 0 2 6*t] * 1/6 * M * p_matrix

            if(with_frames)
                rs[:, a] = [1 t t^2 t^3] * 1/6 * M * spline.minor_ribbon_a[:, i:i+3]'
            end
            a+=1
        end
        a-=1
        i += 1
        
    end



    #####

    sphere_radius = 0.2
    sphere_mesh = discretize(Sphere{3, Float64}((0,0,0), sphere_radius), RegularDiscretization(6))
    debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (200, 200, 200)), eachcol(points)))
    export_mesh_to_ply("main_points.ply", debug_mesh)

    debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (0, 0, 0)), eachcol(points[:, 1:10])))
    export_mesh_to_ply("start.ply", debug_mesh)

    #log_info(misc, "Typ", typeof(points), axes(points))

    if(with_frames)

        debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (100, 0, 0)), eachcol(rs)))
        export_mesh_to_ply("outer_points.ply", debug_mesh)

    #     temp = vcat([itp.(sample_points)' for itp in spline.outer_ribbon_2]...)
    #     debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (0, 100, 0)), eachcol(temp)))
    #     export_mesh_to_ply("outer_points_2.ply", debug_mesh)

        ss = similar(rs)
        normalized_tangents = velocities

        for i=axes(rs, 2)
            normalized_tangents[:, i] /= norm(normalized_tangents[:, i])

            # rs vectors are difference between sampled spline point and the outer spline
            rs[:, i] -= points[:, i]
            #project rs onto plane that is perpendicular to velocities
            rs[:, i] = rs[:, i]/norm(rs[:, i])
            rs[:, i] = rs[:, i] - dot(rs[:, i], normalized_tangents[:, i]) / dot(normalized_tangents[:, i], normalized_tangents[:, i]) * normalized_tangents[:, i]
            rs[:, i] = rs[:, i]/norm(rs[:, i])

            # third axis is perpendicular to the tangent and r
            ss[:, i] = cross(normalized_tangents[:, i], rs[:, i])
        end
        return points, velocities, accelerations, normalized_tangents, rs, ss
    end

    return points, velocities, accelerations
end