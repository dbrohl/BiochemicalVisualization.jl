@enum Part types circle_index_correction time_info damaged_mesh gpu frame_rotation misc point_filter

printed_parts = [gpu, time_info, types, frame_rotation, point_filter, misc] #time_info, circle_index_correction

function log_info(part::Part, args...; separator=" ")
    if(part ∈ printed_parts)
        elems = []
        for a in args
            push!(elems, a)
            push!(elems, separator)
        end
        println(elems...)
    end
end