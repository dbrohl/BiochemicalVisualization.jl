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
    return filter_points(result_points, result_velocities, result_accs)
end

function filter_points(points, velocities, accs)
    min_acc = minimum(norm, eachcol(accs))
    max_acc = maximum(norm, eachcol(accs))

    res_points = Matrix(undef, 3, 0)
    res_vels = Matrix(undef, 3, 0)
    res_accs = Matrix(undef, 3, 0)
    println(min_acc, " ", max_acc)
    for i=axes(points, 2)
        relative_acc = (norm(accs[:, i])-min_acc)/(max_acc-min_acc)
        random_bound = relative_acc*0.8 + 0.2

        r = rand(Float32)

        if(i<10)
            println("relative: $relative_acc, in  [0.5, 1]: $random_bound, r=$r => keep?: $(random_bound>=r)")
        end
        if(i==1 || i==axes(points, 2) || random_bound>=r)
            res_points = [res_points points[:, i]]
            res_vels = [res_vels velocities[:, i]]
            res_accs = [res_accs accs[:, i]]
        end
    end
    println("before: $(size(points, 2)) after: $(size(res_points, 2))")
    return res_points, res_vels, res_accs
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

# TODO find duplicated calculations, granularity for spline
