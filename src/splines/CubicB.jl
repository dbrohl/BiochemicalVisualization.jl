mutable struct CubicB{T <: Real}
    controlPointStrategy
    controlPoints::Matrix{T} # 3 rows, n cols
    minorControlPoints::Union{Matrix{T}, Nothing}

    point_to_residue_indices::Vector{Int}

    num_points_per_resolution::Dict{Int, Vector}
    sample_mapping_per_resolution::Dict{Int, Vector}

    function CubicB(chain::Chain{T}, control_point_strategy) where T #TODO correct first and last points?
        if(control_point_strategy==ControlPoints.C_ALPHA)
            points, point_to_residue_indices = c_alphas_to_points(chain)
            if(length(points)<2)
                throw(ErrorException("too few ($(length(points))) c_alpha atoms to compute spline"))
            end

            first_point = points[1] - (points[2]-points[1])
            last_point = points[end] + (points[end]-points[end-1])
            points = hcat(first_point, points..., last_point)

            prepend!(point_to_residue_indices, point_to_residue_indices[1])
            push!(point_to_residue_indices, point_to_residue_indices[end])

            new{T}(control_point_strategy, points, nothing, point_to_residue_indices, Dict(), Dict())
        elseif(control_point_strategy==ControlPoints.MID_POINTS)
            major_points, minor_points, point_to_residue_indices = generate_points_carson_bugg(chain, true)

            major_points = [major_points[:, 1] major_points major_points[:, end]] #alloc
            minor_points = [minor_points[:, 1] minor_points minor_points[:, end]]

            prepend!(point_to_residue_indices, point_to_residue_indices[1])
            push!(point_to_residue_indices, point_to_residue_indices[end])
            
            new{T}(control_point_strategy, major_points, minor_points, point_to_residue_indices, Dict(), Dict())
        else
            throw(ArgumentError("$control_point_strategy not implemented for CubicBSpline"))
        end
    end
end

function calculate_points(spline::CubicB, resolution)
    dict_key = Int(round(resolution*1000))
    return evaluate_generic_quadruple_spline(spline.controlPoints, num_points(spline, resolution), compute_cubicb_quadruple), sample_to_fragment_index_mapping(spline, resolution)
end

function calculate_velocities(spline::CubicB, resolution)
    return evaluate_generic_quadruple_spline(spline.controlPoints, num_points(spline, resolution), compute_cubicb_quadruple_derivative)
end

function calculate_minor_points(spline::CubicB, resolution)
    return evaluate_generic_quadruple_spline(spline.minorControlPoints, num_points(spline, resolution), compute_cubicb_quadruple)
end

function compute_cubicb_quadruple((P0, P1, P2, P3)::NTuple{4, AbstractVector{T}}, num_points::Int) where T

    result_points = Matrix{T}(undef, 3, num_points)

    M = [1 4 1 0
    -3 0 3 0
    3 -6 3 0
    -1 3 -3 1]
    p_matrix = [P0 P1 P2 P3]'
    fixed_part = T(1/6) * M * p_matrix

    sampling_range = range(T(0), T(1), num_points)

    t_matrix::Matrix{T} = [1 0 0 0]
    for (i, t) in enumerate(sampling_range)
        t_matrix[1, 2:4] = [t t^2 t^3]
        result_points[:, i] = t_matrix * fixed_part
    end
    return result_points
end

function compute_cubicb_quadruple_derivative((P0, P1, P2, P3)::NTuple{4, AbstractVector{T}}, num_points::Int) where T
    result_points = Matrix{T}(undef, 3, num_points)

    M = [1 4 1 0
    -3 0 3 0
    3 -6 3 0
    -1 3 -3 1]
    p_matrix = [P0 P1 P2 P3]'
    fixed_part = T(1/6) * M * p_matrix

    sampling_range = range(T(0), T(1), num_points)
    t_matrix::Matrix{T} = [0 1 0 0]
    for (i, t) in enumerate(sampling_range)
        t_matrix[1, 3:4] = [2*t 3*t^2]
        result_points[:, i] = t_matrix * fixed_part
    end
    return result_points
end