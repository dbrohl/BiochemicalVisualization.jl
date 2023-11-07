mutable struct CatmullRom
    controlPointStrategy
    controlPoints::Matrix # 3 rows, n cols
    minorControlPoints::Union{Matrix, Nothing}

    point_to_residue_indices::Vector

    num_points_per_resolution::Dict{Int, Vector}
    sample_mapping_per_resolution::Dict{Int, Vector}

    function CatmullRom(chain::BiochemicalAlgorithms.Chain, control_point_strategy)
        if(control_point_strategy==ControlPoints.C_ALPHA)
            points, point_to_residue_indices = c_alphas_to_points(chain)
            if(length(points)<2)
                throw(ErrorException("too few ($(length(points))) c_alpha atoms to compute spline"))
            end

            # add first and last dummy point
            # control points cannot be the same (otherwise the sampling produces NaN values)
            first_point = points[1] - (points[2]-points[1])
            last_point = points[end] + (points[end]-points[end-1])
            points = hcat(first_point, points..., last_point)

            prepend!(point_to_residue_indices, point_to_residue_indices[1])
            push!(point_to_residue_indices, point_to_residue_indices[end])

            new(control_point_strategy, points, nothing, point_to_residue_indices, Dict(), Dict())
        elseif(control_point_strategy==ControlPoints.MID_POINTS)
            major_points, minor_points, point_to_residue_indices = generate_points_carson_bugg(chain, false)

            # add first and last dummy point
            # control points cannot be the same (otherwise the sampling produces NaN values)
            major_points = [major_points[:, 1]-(major_points[:,2]-major_points[:,1]) major_points major_points[:, end]-(major_points[:, end]-major_points[:, end-1])]
            minor_points = [minor_points[:, 1]-(minor_points[:,2]-minor_points[:,1]) minor_points minor_points[:, end]-(minor_points[:, end]-minor_points[:, end-1])]

            prepend!(point_to_residue_indices, point_to_residue_indices[1])
            push!(point_to_residue_indices, point_to_residue_indices[end])
            
            new(control_point_strategy, major_points, minor_points, point_to_residue_indices, Dict(), Dict())
        else
            throw(ArgumentError("$control_point_strategy not implemented for CatmullRomSpline"))
        end
    end
end

function calculate_points(spline::CatmullRom, resolution)
    dict_key = Int(round(resolution*1000))
    return evaluate_generic_quadruple_spline(spline.controlPoints, num_points(spline, resolution), compute_catmull_rom_quadruple), sample_to_fragment_index_mapping(spline, resolution)
end

function calculate_velocities(spline::CatmullRom, resolution)
    return evaluate_generic_quadruple_spline(spline.controlPoints, num_points(spline, resolution), compute_catmull_rom_quadruple_derivative)
end

function calculate_minor_points(spline::CatmullRom, resolution)
    return evaluate_generic_quadruple_spline(spline.minorControlPoints, num_points(spline, resolution), compute_catmull_rom_quadruple)
end



# Code adapted from https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline#Code_example_in_Python (Last access: 24.07.2023)
function compute_catmull_rom_quadruple((P0, P1, P2, P3), num_points)
    t0 = 0
    t1 = tRecursion(P1, P0, t0)
    t2 = tRecursion(P2, P1, t1)
    t3 = tRecursion(P3, P2, t2)

    result_points = Matrix(undef, 3, num_points)

    ts = collect(range(t1, t2, num_points))

    for (i, t) in enumerate(ts)
        A1 = @. (t1-t)/(t1-t0) * P0 + (t-t0)/(t1-t0) * P1
        A2 = @. (t2-t)/(t2-t1) * P1 + (t-t1)/(t2-t1) * P2
        A3 = @. (t3-t)/(t3-t2) * P2 + (t-t2)/(t3-t2) * P3

        B1 = @. (t2-t)/(t2-t0) * A1 + (t-t0)/(t2-t0) * A2
        B2 = @. (t3-t)/(t3-t1) * A2 + (t-t1)/(t3-t1) * A3

        C = @. (t2-t)/(t2-t1) * B1 + (t-t1)/(t2-t1) * B2

        result_points[:, i] = C
    end
    return result_points
end

function compute_catmull_rom_quadruple_derivative((P0, P1, P2, P3), num_points)
    t0 = 0
    t1 = tRecursion(P1, P0, t0)
    t2 = tRecursion(P2, P1, t1)
    t3 = tRecursion(P3, P2, t2)

    result_velocities = Matrix(undef, 3, num_points)
    ts = collect(range(t1, t2, num_points))

    for (i, t) in enumerate(ts)
        A1 = @. (t1-t)/(t1-t0) * P0 + (t-t0)/(t1-t0) * P1
        A2 = @. (t2-t)/(t2-t1) * P1 + (t-t1)/(t2-t1) * P2
        A3 = @. (t3-t)/(t3-t2) * P2 + (t-t2)/(t3-t2) * P3

        B1 = @. (t2-t)/(t2-t0) * A1 + (t-t0)/(t2-t0) * A2
        B2 = @. (t3-t)/(t3-t1) * A2 + (t-t1)/(t3-t1) * A3

        # ----- first derivative -----
        A1v = @. -1/(t1-t0)*P0 + 1/(t1-t0)*P1
        A2v = @. -1/(t2-t1)*P1 + 1/(t2-t1)*P2
        A3v = @. -1/(t3-t2)*P2 + 1/(t3-t2)*P3

        B1v = @. -1/(t2-t0)*A1 + (t2-t)/(t2-t0)*A1v + 1/(t2-t0)*A2 + (t-t0)/(t2-t0)*A2v
        B2v = @. -1/(t3-t1)*A2 + (t3-t)/(t3-t1)*A2v + 1/(t3-t1)*A3 + (t-t1)/(t3-t1)*A3v
        Cv  = @. -1/(t2-t1)*B1 + (t2-t)/(t2-t1)*B1v + 1/(t2-t1)*B2 + (t-t1)/(t2-t1)*B2v

        result_velocities[:, i] = Cv
    end
    return result_velocities
end

function tRecursion(pCurr, pPrev, tPrev)

    distance = norm(pCurr .- pPrev)
    return distance^0.5 + tPrev
end
