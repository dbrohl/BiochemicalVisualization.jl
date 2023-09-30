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
function filter_points_threshold(points, tangents)
    remaining_indices = []
    for i=axes(points, 2)
        if(i==firstindex(points, 2) || i==lastindex(points, 2))
            println("Found start or end ", i)
            push!(remaining_indices, i)
            continue
        end


        angle_between_frames = acos(dot(tangents[:, remaining_indices[end]], tangents[:, i])) #green: acos(dot(normals[:, remaining_indices[end]], normals[:, i]))
        #log_info(point_filter, i, angle_between_frames/(2*π)*360, normals[:, remaining_indices[end]], normals[:, i], norm(normals[:, remaining_indices[end]]), norm(normals[:, i]))
        if(abs(angle_between_frames)> 5/360*2*π)
            push!(remaining_indices, i)
        end
    end

    log_info(point_filter, "before: $(size(points, 2)) after: $(length(remaining_indices))")
    return remaining_indices
end