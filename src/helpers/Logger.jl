@enum Part types circle_index_correction time_info damaged_mesh gpu frame_rotation

printed_parts = [gpu, time_info, types, frame_rotation] #time_info, circle_index_correction

function log_info(part::Part, args...)
    if(part ∈ printed_parts)
        println(args...)
    end
end