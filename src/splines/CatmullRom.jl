export CatmullRom

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

    result = Matrix(undef, 3, sum(num_points)-length(num_points)+1)

    i = 1
    a = 1
    while i+3 <= size(spline.controlPoints, 2) # loop over quadruple of controlPoints
        
        points = compute_quadruple((spline.controlPoints[:, i], spline.controlPoints[:, i+1], spline.controlPoints[:, i+2], spline.controlPoints[:, i+3]), num_points[i])
        result[:, a:a+num_points[i]-1] = points

        a += num_points[i]-1
        i += 1
        
    end
    return result
end

function compute_quadruple((P0, P1, P2, P3), num_points)
    t0 = 0
    t1 = tRecursion(P1, P0, t0)
    t2 = tRecursion(P2, P1, t1)
    t3 = tRecursion(P3, P2, t2)

    result = Matrix(undef, 3, num_points)
    ts = collect(range(t1, t2, num_points))

    for (i, t) in enumerate(ts)
        A1 = @. (t1-t)/(t1-t0) * P0 + (t-t0)/(t1-t0) * P1
        A2 = @. (t2-t)/(t2-t1) * P1 + (t-t1)/(t2-t1) * P2
        A3 = @. (t3-t)/(t3-t2) * P2 + (t-t2)/(t3-t2) * P3

        B1 = @. (t2-t)/(t2-t0) * A1 + (t-t0)/(t2-t0) * A2
        B2 = @. (t3-t)/(t3-t1) * A2 + (t-t1)/(t3-t1) * A3

        C = @. (t2-t)/(t2-t1) * B1 + (t-t1)/(t2-t1) * B2
        result[:, i] = C
    end
    return result
end

function tRecursion(pCurr, pPrev, tPrev)

    distance = norm(pCurr .- pPrev)
    return distance^0.5 + tPrev
end

# TODO find duplicated calculations, granularity for spline
