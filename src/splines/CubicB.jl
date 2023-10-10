mutable struct CubicB
    controlPoints::AbstractArray
    interpolations::Union{Vector{SplineInterpolation}, Nothing}

    function CubicB(chain)
        n_ribbons = 3
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
        # sphere_radius = 0.2
        # sphere_mesh = discretize(Sphere{3, Float64}((0,0,0), sphere_radius), RegularDiscretization(6))
        # debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a.r)...)((sphere_mesh)), (0, 0, 0)), c_alphas))
        # export_mesh_to_ply("c_alphas.ply", debug_mesh)
        # debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a.r)...)((sphere_mesh)), (0, 0, 255)), oxygens))
        # export_mesh_to_ply("oxygens.ply", debug_mesh)

        controlPoints = Array{Float64, 3}(undef, n_ribbons, 3, length(c_alphas)-1)
        for i=1:length(c_alphas)-1
            A = c_alphas[i+1].r - c_alphas[i].r
            B = oxygens[i].r - c_alphas[i].r
            C = cross(A, B)
            D = cross(C, A)
    
            C = C/norm(C)
            D = D/norm(D)
    
            P = c_alphas[i].r + 0.5 * A
            if(structures[i]==SecondaryStructure.HELIX || structures[i+1]==SecondaryStructure.HELIX) # TODO correct condition?
                # TODO translate
            end
            D *= 0.5 * ribbon_width[structures[i]]
    
            for (j, t) in enumerate(range(-1, 1, n_ribbons))
                controlPoints[j, :, i] = P + t*D
            end
    
        end

        # debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (255, 0, 0)), eachcol(controlPoints[2, :, :])))
        # export_mesh_to_ply("control_points.ply", debug_mesh)

        new(controlPoints, nothing)
    end
end

function (spline::CubicB)(resolution)
    
    sample_points = []
    for i=1:size(spline.controlPoints, 3)-1
        distance = norm(spline.controlPoints[2, :, i] .- spline.controlPoints[2, :, i+1])
        num_points = max(2, convert(Int, ceil(resolution * 1 * distance)))

        if(i==1)
            push!(sample_points, range(i, i+1, num_points)...)
        else
            push!(sample_points, range(i, i+1, num_points)[2:end]...)
        end
    end

    if(spline.interpolations === nothing)
        knots = 1:size(spline.controlPoints, 3)
        spline.interpolations = [interpolate(knots, spline.controlPoints[2, 1, :], BSplineOrder(4)),
                                 interpolate(knots, spline.controlPoints[2, 2, :], BSplineOrder(4)),
                                 interpolate(knots, spline.controlPoints[2, 3, :], BSplineOrder(4))]
    end

    points = vcat([itp.(sample_points)' for itp in spline.interpolations]...)
    velocities = vcat([(Derivative(1)*itp).(sample_points)' for itp in spline.interpolations]...)
    accelerations = vcat([(Derivative(2)*itp).(sample_points)' for itp in spline.interpolations]...)

    # sphere_radius = 0.2
    # sphere_mesh = discretize(Sphere{3, Float64}((0,0,0), sphere_radius), RegularDiscretization(6))
    # debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (0, 255, 0)), eachcol(points)))
    # export_mesh_to_ply("points.ply", debug_mesh)

    #log_info(misc, "Typ", typeof(points), axes(points))
    return points, velocities, accelerations
end
