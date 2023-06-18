export CatmullRom

struct CatmullRom
    controlPoints::AbstractArray{Vec} # vectors support +,*,... (Points do not)

    function CatmullRom(points::AbstractVector{T}) where T<:Point # add first and last point
        points = map(p -> p.coords, points)
        CatmullRom(points)
    end

    function CatmullRom(points::AbstractVector{T}) where T<:Vec # add first and last point
        prepend!(points, [points[1]-(points[2]-points[1])])
        push!(points, points[end]+(points[end]-points[end-1]))
        new(points)
    end
end

function (spline::CatmullRom)() #t::Real
    result = []
    i=1
    while i+3 <= length(spline.controlPoints)
        push!(result, 
            compute_quadruple(spline.controlPoints[i], spline.controlPoints[i+1], spline.controlPoints[i+2], spline.controlPoints[i+3])...)
        i+=1
    end
    return result
end

function compute_quadruple(P0, P1, P2, P3)
    t0 = 0
    t1 = tRecursion(P1, P0, t0)
    t2 = tRecursion(P2, P1, t1)
    t3 = tRecursion(P3, P2, t2)

    result = []
    ts = collect(range(t1, t2, 100))

    for t in ts
        A1 = (t1-t)/(t1-t0) * P0 + (t-t0)/(t1-t0) * P1
        A2 = (t2-t)/(t2-t1) * P1 + (t-t1)/(t2-t1) * P2
        A3 = (t3-t)/(t3-t2) * P2 + (t-t2)/(t3-t2) * P3

        B1 = (t2-t)/(t2-t0) * A1 + (t-t0)/(t2-t0) * A2
        B2 = (t3-t)/(t3-t1) * A2 + (t-t1)/(t3-t1) * A3

        C = (t2-t)/(t2-t1) * B1 + (t-t1)/(t2-t1) * B2
        push!(result, C)
    end
    return result
end

function tRecursion(pCurr, pPrev, tPrev)

    distance = sqrt(sum((pCurr .- pPrev) .^ 2))
    return distance^0.5 + tPrev
end

# TODO credit, find duplicated calculations, granularity for spline
