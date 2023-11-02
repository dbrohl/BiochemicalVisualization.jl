# Multiple methods to filter the sampled points of a spline

function no_filter(points)
    return axes(points, 2)
end

# Smaller lengths of the acceleration vectors correspond to a higher probability that the point gets deleted. 
function filter_points_stoch(points, accs)
    min_acc = minimum(norm, eachcol(accs))
    max_acc = maximum(norm, eachcol(accs))

    remaining_indices = []
    for i=axes(points, 2)
        relative_acc = (norm(accs[:, i])-min_acc)/(max_acc-min_acc)
        random_bound = relative_acc*0.8 + 0.2

        r = rand(Float32)

        if(i==firstindex(points, 2) || i==lastindex(points, 2) || random_bound>=r)
            push!(remaining_indices, i)
        end
    end
    log_info(point_filter, "before: $(size(points, 2)) after: $(length(remaining_indices))")
    return remaining_indices
end

# Whenever the angle between the last selected and the current tangent is too large, the current circle is added. 
# Tangents should be normalized!
# function filter_points_threshold(points, tangents)
#     remaining_indices = []

#     for i=axes(points, 2)
#         if(i==firstindex(points, 2) || i==lastindex(points, 2))
#             push!(remaining_indices, i)
#             continue
#         end

#         dot_prod = dot(tangents[:, remaining_indices[end]], tangents[:, i])
#         if(abs(dot_prod)>1 || abs(acos(dot_prod))> 5/360*2*π) # Numerical issues could lead to a DomainError when the dot_product is slightly larger than 1. 
#             #log_info(point_filter, i, acos(dot_prod))
#             push!(remaining_indices, i)
#         end
#     end

#     log_info(point_filter, "before: $(size(points, 2)) after: $(length(remaining_indices))")
#     return remaining_indices
# end

# assumes at least 1 in fixed indices
function filter_points_threshold(points, q, r, s, fixed_indices, colors=nothing)
    remaining_indices = []
    a = 0
    b = 0
    c = 0

    degree_threshold = 10
    color_degree_threshold = 45

    for i=axes(points, 2)
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