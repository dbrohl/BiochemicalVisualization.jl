struct CatmullRom
    controlPoints::AbstractMatrix # 3 rows, n cols

    function CatmullRom(points::AbstractMatrix)
        points = hcat(points[:, 1] .- (points[:, 2] .- points[:, 1]), points, points[:, end] .+ (points[:, end] .- points[:, end-1]))
        new(points)
    end
end

# Code adapted from https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline#Code_example_in_Python (Last access: 24.07.2023)
function (spline::CatmullRom)(resolution)
    num_points = []
    i = 1
    while i+3 <= size(spline.controlPoints, 2)
        distance = norm(spline.controlPoints[:, i+1] .- spline.controlPoints[:, i+2])
        push!(num_points, max(2, convert(Int, ceil(resolution * 2 * distance))))
        i += 1
    end

    result_points = Matrix(undef, 3, sum(num_points)-length(num_points)+1)
    result_velocities = Matrix(undef, 3, sum(num_points)-length(num_points)+1)
    result_accs = Matrix(undef, 3, sum(num_points)-length(num_points)+1)

    i = 1
    a = 1
    while i+3 <= size(spline.controlPoints, 2) # loop over quadruple of controlPoints
        
        points, velocities, accs = compute_quadruple((spline.controlPoints[:, i], spline.controlPoints[:, i+1], spline.controlPoints[:, i+2], spline.controlPoints[:, i+3]), num_points[i])
        result_points[:, a:a+num_points[i]-1] = points
        result_velocities[:, a:a+num_points[i]-1] = velocities
        result_accs[:, a:a+num_points[i]-1] = accs
        a += num_points[i]-1
        i += 1
        
    end

    initial_tangent = result_velocities[:, 1]./norm(result_velocities[:, 1])
    if(approx_zero(initial_tangent[2]) && approx_zero(initial_tangent[3]))
        temp = [1; 0; 0]
    else
        temp = [0; 1; 0]
    end
    initial_r = cross(initial_tangent, temp)
    initial_r = initial_r ./ norm(initial_r)
    ts, rs, ss = rmf(result_points, result_velocities, initial_r, cross(initial_tangent, initial_r))

    return result_points, result_velocities, result_accs, (ts, rs, ss)
end



# based on W. Wang, B. Jüttler, D. Zheng, and Y. Liu, ‘Computation of rotation minimizing frames’, doi: 10.1145/1330511.1330513.
function rmf(points, tangents, r0, s0)

    ts = Matrix(undef, 3, size(points, 2))
    rs = Matrix(undef, 3, size(points, 2))
    ss = Matrix(undef, 3, size(points, 2))

    for (i, col) in enumerate(eachcol(tangents))
        ts[:, i] = col ./ norm(col)
    end

    rs[:, 1] = r0
    ss[:, 1] = s0

    for i=1:size(points, 2)-1
        v1 = points[:, i+1] - points[:, i]
        c1 = dot(v1,v1)
        r_i_L = rs[:, i] - (2/c1) * dot(v1, rs[:, i]) * v1
        t_i_L = ts[:, i] - (2/c1) * dot(v1, ts[:, i]) * v1
        v2 = ts[:, i+1] - t_i_L
        c2 = dot(v2, v2)
        rs[:, i+1] = r_i_L - (2/c2) * dot(v2, r_i_L) * v2
        rs[:, i+1] = rs[:, i+1] ./ norm(rs[:, i+1])
        ss[:, i+1] = cross(ts[:, i+1], rs[:, i+1])
    end
    return ts, rs, ss

end



function compute_quadruple((P0, P1, P2, P3), num_points)
    t0 = 0
    t1 = tRecursion(P1, P0, t0)
    t2 = tRecursion(P2, P1, t1)
    t3 = tRecursion(P3, P2, t2)

    result_points = Matrix(undef, 3, num_points)
    result_velocities = Matrix(undef, 3, num_points)
    result_accs = Matrix(undef, 3, num_points)
    ts = collect(range(t1, t2, num_points))

    for (i, t) in enumerate(ts)
        A1 = @. (t1-t)/(t1-t0) * P0 + (t-t0)/(t1-t0) * P1
        A2 = @. (t2-t)/(t2-t1) * P1 + (t-t1)/(t2-t1) * P2
        A3 = @. (t3-t)/(t3-t2) * P2 + (t-t2)/(t3-t2) * P3

        B1 = @. (t2-t)/(t2-t0) * A1 + (t-t0)/(t2-t0) * A2
        B2 = @. (t3-t)/(t3-t1) * A2 + (t-t1)/(t3-t1) * A3

        C = @. (t2-t)/(t2-t1) * B1 + (t-t1)/(t2-t1) * B2

        # ----- first derivative -----
        A1v = @. -1/(t1-t0)*P0 + 1/(t1-t0)*P1
        A2v = @. -1/(t2-t1)*P1 + 1/(t2-t1)*P2
        A3v = @. -1/(t3-t2)*P2 + 1/(t3-t2)*P3

        B1v = @. -1/(t2-t0)*A1 + (t2-t)/(t2-t0)*A1v + 1/(t2-t0)*A2 + (t-t0)/(t2-t0)*A2v
        B2v = @. -1/(t3-t1)*A2 + (t3-t)/(t3-t1)*A2v + 1/(t3-t1)*A3 + (t-t1)/(t3-t1)*A3v
        Cv  = @. -1/(t2-t1)*B1 + (t2-t)/(t2-t1)*B1v + 1/(t2-t1)*B2 + (t-t1)/(t2-t1)*B2v


        # ----- second derivative -----
        A1a = 0
        A2a = 0
        A3a = 0

        B1a = @. -1/(t2-t0)*A1v + -1/(t2-t0)*A1v + (t2-t)/(t2-t0)*A1a + 1/(t2-t0)*A2v + 1/(t2-t0)*A2v + (t-t0)/(t2-t0)*A2a
        B2a = @. -1/(t3-t1)*A2v + -1/(t3-t1)*A2v + (t3-t)/(t3-t1)*A2a + 1/(t3-t1)*A3v + 1/(t3-t1)*A3v + (t-t1)/(t3-t1)*A3a
        Ca  = @. -1/(t2-t1)*B1v + -1/(t2-t1)*B1v + (t2-t)/(t2-t1)*B1a + 1/(t2-t1)*B2v + 1/(t2-t1)*B2v + (t-t1)/(t2-t1)*B2a


        # ----- third derivative -----
        A1j = 0
        A2j = 0
        A3j = 0

        B1j = @. -1/(t2-t0)*A1a + -1/(t2-t0)*A1a + -1/(t2-t0)*A1a + (t2-t)/(t2-t0)*A1j + 1/(t2-t0)*A2a + 1/(t2-t0)*A2a + 1/(t2-t0)*A2a + (t-t0)/(t2-t0)*A2j
        B2j = @. -1/(t3-t1)*A2a + -1/(t3-t1)*A2a + -1/(t3-t1)*A2a + (t3-t)/(t3-t1)*A2j + 1/(t3-t1)*A3a + 1/(t3-t1)*A3a + 1/(t3-t1)*A3a + (t-t1)/(t3-t1)*A3j
        Cj  = @. -1/(t2-t1)*B1a + -1/(t2-t1)*B1a + (-1)/(t2-t1)*B1a + (t2-t)/(t2-t1)*B1j + 1/(t2-t1)*B2a + 1/(t2-t1)*B2a + 1/(t2-t1)*B2a + (t-t1)/(t2-t1)*B2j
        

        result_points[:, i] = C
        result_velocities[:, i] = Cv
        result_accs[:, i] = Ca
    end
    return result_points, result_velocities, result_accs
end

function tRecursion(pCurr, pPrev, tPrev)

    distance = norm(pCurr .- pPrev)
    return distance^0.5 + tPrev
end
