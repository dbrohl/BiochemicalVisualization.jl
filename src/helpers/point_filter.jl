# Multiple methods to filter the sampled points of a spline

function no_filter(points)
    return axes(points, 2)
end

function filter_points_stoch(points, accs)
    min_acc = minimum(norm, eachcol(accs))
    max_acc = maximum(norm, eachcol(accs))
    println(min_acc, " ", max_acc)

    remaining_indices = []
    for i=axes(points, 2)
        relative_acc = (norm(accs[:, i])-min_acc)/(max_acc-min_acc)
        random_bound = relative_acc*0.8 + 0.2

        r = rand(Float32)

        if(i==1 || i==axes(points, 2) || random_bound>=r)
            push!(remaining_indices, i)
        end
    end
    println("before: $(size(points, 2)) after: $(len(remaining_indices))")
    return remaining_indices
end