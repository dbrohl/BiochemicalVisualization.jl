"""
Returns the indices of points that should not be removed. 
Whenever the angle between the last selected and the current tangent is too large, the current index is added. 

* When colors!=nothing, too large hue distances are prevented as well.
* Vectors in q and r should be normalized!
* fixed_indices contains all indices that cannot be removed and will definitely be contained in teh return value. 
"""
function filter_points_threshold(q::Matrix{T}, r::Matrix{T}, fixed_indices::Vector{Int}, colors::Union{Nothing, Vector{NTuple{3, Int}}}=nothing) where T
    if(length(fixed_indices)==0)
        fixed_indices = [1]
    end
    remaining_indices = []
    a = 0
    b = 0
    c = 0

    degree_threshold = 5
    color_degree_threshold = 45

    for i=axes(q, 2)
        if(i ∈ fixed_indices)
            push!(remaining_indices, i)
        else
            dot_prod_q = dot(q[:, remaining_indices[end]], q[:, i])
            if(abs(dot_prod_q)>1 || abs(acos(dot_prod_q))> degree_threshold/360*2*π) # Numerical issues could lead to a DomainError when the dot_product is slightly larger than 1.
                a+=1 
            end

            dot_prod_r = dot(r[:, remaining_indices[end]], r[:, i])
            if(abs(dot_prod_r)>1 || abs(acos(dot_prod_r))> degree_threshold/360*2*π) # Numerical issues could lead to a DomainError when the dot_product is slightly larger than 1.
                b+=1 
            end

            large_color_distance = false
            if(colors!==nothing)
                colA = convert(HSL, RGB(colors[remaining_indices[end]]./255...))
                colB = convert(HSL, RGB(colors[i]./255...))

                if(abs(colA.h-colB.h)>color_degree_threshold) 
                    large_color_distance = true
                    c += 1
                end
            end

            if(abs(dot_prod_q)>1 || abs(acos(dot_prod_q))> degree_threshold/360*2*π 
                || abs(dot_prod_r)>1 || abs(acos(dot_prod_r))> degree_threshold/360*2*π
                || large_color_distance)
                push!(remaining_indices, i)
            end

        end
    end
    log_info(point_filter, "filtered by tangents: $a, filtered by normals: $b, filtered by color: $c")
    return remaining_indices
end