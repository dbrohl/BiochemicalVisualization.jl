@enum Part types circle_index_correction time_info damaged_mesh gpu frame_rotation misc point_filter config

printed_parts = []#gpu, time_info, types, point_filter, misc, frame_rotation]

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

function log_warning(part::Part, args...; separator=" ")
    elems = []
    for a in args
        push!(elems, a)
        push!(elems, separator)
    end
    println(elems...)
end