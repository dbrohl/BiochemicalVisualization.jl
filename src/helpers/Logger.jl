@enum Part types circle_index_correction time_info damaged_mesh gpu frame_rotation misc point_filter config extra_frames short_chains

printed_parts = [time_info, types]#gpu, time_info, types, point_filter, misc, frame_rotation]

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

function log_sizes(part::Part, arrays...)
    if(part ∈ printed_parts)
        elems = []
        for a in arrays
            if(a===nothing)
                push!(elems, "nothing")
            else
                push!(elems, size(a))
            end
            push!(elems, " ")
        end
        println(elems...)
    end
end

function log_warning(args...; separator=" ")
    elems = []
    for a in args
        push!(elems, a)
        push!(elems, separator)
    end
    println(elems...)
end