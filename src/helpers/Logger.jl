@enum Part types circle_index_correction time_info damaged_mesh gpu

printed_parts = [gpu, time_info, types] #time_info, circle_index_correction

function log_info(part::Part, args...)
    if(part ∈ printed_parts)
        println(args...)
    end
end