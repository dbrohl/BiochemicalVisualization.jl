"""
Returns the indices of points that should not be removed. 
Whenever the angle between the last selected and the current tangent is too large, the current index is added. 

* When colors!=nothing, too large hue distances are prevented as well.
* Vectors in q and r should be normalized!
* fixed_indices contains all indices that cannot be removed and will definitely be contained in the return value. 
"""
function filter_points_threshold(q::Matrix{T}, r::Matrix{T}, fixed_indices::AbstractVector{Int}, colors::Union{Nothing, AbstractVector{NTuple{3, Int}}}=nothing) where T
    remaining_indices = Set(fixed_indices)
    push!(remaining_indices, 1)
    last_remaining_index = 1

    degree_threshold = 5
    radian_threshold = degree_threshold/360*2*π 
    color_degree_threshold = 45

    dot_prod_q = T(0)
    dot_prod_r = T(0)
    large_color_distance = false

    for i=axes(q, 2)
        if(i ∈ remaining_indices)
            last_remaining_index = i
        else
            dot_prod_q = @views dot(q[:, last_remaining_index], q[:, i]) # alloc
            dot_prod_r = @views dot(r[:, last_remaining_index], r[:, i])

            large_color_distance = false
            if(colors!==nothing)
                colA = convert(HSL, RGB(colors[last_remaining_index]./255...)).h
                colB = convert(HSL, RGB(colors[i]./255...)).h

                large_color_distance = abs(colA-colB)>color_degree_threshold
            end

            if(abs(dot_prod_q)>1 || abs(acos(dot_prod_q))> radian_threshold
                || abs(dot_prod_r)>1 || abs(acos(dot_prod_r))> radian_threshold
                || large_color_distance)
                push!(remaining_indices, i)
                last_remaining_index = i
            end

        end
    end
    return remaining_indices
end