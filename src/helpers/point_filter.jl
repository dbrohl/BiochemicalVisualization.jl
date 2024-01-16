"""
Returns the indices and the number of points that should not be removed. 
Whenever the angle between the last selected and the current tangent is too large, the current index is added. 

* If with_color, a linear hsv interpolation is assumed. 
  When colors change too much (in this interpolation), the current index is added to the result. 
* Vectors in q and r should be normalized!
* fixed_indices contains all indices that cannot be removed and will definitely be contained in the return value. 
"""
function filter_points_threshold(q::Matrix{T}, r::Matrix{T}, fixed_indices::AbstractVector{Int}; with_color::Bool=false) where T

    target_indices = fill(-1, size(q, 2))
    sort!(fixed_indices)

    degree_threshold = 5
    radian_threshold = degree_threshold/360*2*π 
    color_threshold = 0.125 # max part of the rainbow that will be interpolated (vs an additional frame to ensure proper colors)

    dot_prod_q = T(0)
    dot_prod_r = T(0)
    large_color_distance = false

    
    last_remaining_index = 1
    a = 1
    for i=axes(q, 2)
        if(i==1 || insorted(i, fixed_indices))
            last_remaining_index = i
            target_indices[i] = a
            a += 1
        else
            dot_prod_q = @views dot(q[:, last_remaining_index], q[:, i]) # alloc
            dot_prod_r = @views dot(r[:, last_remaining_index], r[:, i])

            large_color_distance = with_color && (i-last_remaining_index)/size(q, 2)>color_threshold

            if(abs(dot_prod_q)>1 || abs(acos(dot_prod_q))> radian_threshold
                || abs(dot_prod_r)>1 || abs(acos(dot_prod_r))> radian_threshold
                || large_color_distance)
                last_remaining_index = i
                target_indices[i] = a
                a += 1
            end

        end
    end
    return target_indices, a-1
end