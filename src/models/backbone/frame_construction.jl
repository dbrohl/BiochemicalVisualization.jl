"""
Creates orthonormal frames along a series of points and their tangents. The rotation around the tangent-axis is minimized. 
Returns normalized (tangents, normals, binormals). 

based on W. Wang, B. Jüttler, D. Zheng, and Y. Liu, "Computation of rotation minimizing frames", doi: 10.1145/1330511.1330513.
"""
function rmf(points::Matrix{T}, tangents::Matrix{T}) where T
    ts = Matrix{T}(undef, 3, size(points, 2))
    rs = Matrix{T}(undef, 3, size(points, 2))
    ss = Matrix{T}(undef, 3, size(points, 2))

    for (i, col) in enumerate(eachcol(tangents)) #alloc
        if(approx_zero(norm(col)))
            log_warning("zero length tangent in rmf (index $i/$(size(tangents, 2))): $col")
        end
        @. ts[:, i] = col / norm(col)
    end

    if(approx_zero(ts[2, 1]) && approx_zero(ts[3, 1]))
        temp = [0; 1; 0]
    else
        temp = [1; 0; 0]
    end
    rs[:, 1] = @views cross( ts[:, 1], temp)
    @views rs[:, 1] ./= norm(rs[:, 1])

    ss[:, 1] = cross(ts[:, 1], rs[:, 1])

    @views for i=1:size(points, 2)-1
        v1 = @. points[:, i+1] - points[:, i]
        c1 = dot(v1,v1)
        r_i_L = rs[:, i] .- ((2/c1) * dot(v1, rs[:, i])) .* v1
        t_i_L = ts[:, i] .- ((2/c1) * dot(v1, ts[:, i])) .* v1
        v2 = ts[:, i+1] - t_i_L
        c2 = dot(v2, v2)
        rs[:, i+1] .= r_i_L .- ((2/c2) * dot(v2, r_i_L)) .* v2
        if approx_zero(norm(rs[:, i+1]))
            log_warning("zero length r in rmf (index $(i+1)): $(rs[:, i+1])) $(norm(rs[:, i+1]))")
        end
        rs[:, i+1] ./= norm(rs[:, i+1])
        ss[:, i+1] = cross(ts[:, i+1], rs[:, i+1])
    end
    return ts, rs, ss
end

"""
Creates orthonormal frames along a series of points (derived from the major spline). 
The normals point towards the respective points on the minor spline. 
Returns normalized (tangents, normals, binormals). 

Approch by M. Carson and C. E. Bugg, "Algorithm for ribbon models of proteins", doi: 10.1016/0263-7855(86)80010-8.
"""
function frames_from_two_splines(major_spline_points::Matrix{T}, major_spline_tangents::Matrix{T}, minor_spline_points::Matrix{T}) where T
    ts = similar(major_spline_tangents)
    rs = similar(major_spline_tangents)
    ss = similar(major_spline_tangents)


    @views for i=axes(major_spline_tangents, 2)
        ts[:, i] = major_spline_tangents[:, i] ./ norm(major_spline_tangents[:, i])

        # rs vectors are difference between sampled spline point and the outer spline
        rs[:, i] = minor_spline_points[:, i] .- major_spline_points[:, i]
        #project r onto plane that is perpendicular to tangent
        rs[:, i] ./= norm(rs[:, i])
        rs[:, i] .-= (dot(rs[:, i], ts[:, i]) / dot(ts[:, i], ts[:, i]) .* ts[:, i])
        rs[:, i] ./= norm(rs[:, i])

        # third axis is perpendicular to the tangent and r
        ss[:, i] = cross(ts[:, i], rs[:, i]) #alloc why?
    end
    return ts, rs, ss
end