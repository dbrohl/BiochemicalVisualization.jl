export prepare_backbone_model, generate_chain_mesh, generate_geometry_at_point, benchmark_method # TODO remove

function benchmark_method(ac::System{T}, config::BackboneConfig) where T
    mesh = prepare_backbone_model(ac, config)
    representation = Representation(mesh) 
    return representation
end

function insert_sorted!(array, elem)
    index = searchsortedfirst(array, elem)
    insert!(array, index, elem)
end

function insert_frame(index, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = nothing) # TODO types?
    if(index<1 || index>size(spline_points, 2))
        throw(ArgumentError("index has to be in [1, $(size(spline_points, 2))], but was $index"))
    end
    spline_points = @views hcat(spline_points[:, 1:index], spline_points[:, index], spline_points[:, index+1:end])
    q = @views hcat(q[:, 1:index], q[:, index], q[:, index+1:end])
    r = @views hcat(r[:, 1:index], r[:, index], r[:, index+1:end])
    s = @views hcat(s[:, 1:index], s[:, index], s[:, index+1:end])
    rectangle_widths = @views vcat(rectangle_widths[1:index], rectangle_widths[index], rectangle_widths[index+1:end])
    sample_to_residue_indices = @views vcat(sample_to_residue_indices[1:index], sample_to_residue_indices[index], sample_to_residue_indices[index+1:end])
    if(rainbow_colors !== nothing)
        rainbow_colors = @views vcat(rainbow_colors[1:index], rainbow_colors[index], rainbow_colors[index+1:end])
    end

    return spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors
end

#assumes sorted index_array!
function adjust_indices!(index_array, threshold)
    for i=eachindex(index_array)
        if(index_array[i]>=threshold)
            index_array[i] += 1
        end
    end
end

function compute_frame_widths(fragment_list::Vector{Fragment{T}}, sample_to_residue_indices) where T
    rectangle_widths::Vector{T} = []
    arrow_starts::Vector{Int} = []
    arrow_frame_indices::Vector{Int} = []

    filtered_fragment_list = filter(f -> is_amino_acid(f), fragment_list)
    # determine orientation of looping (arrows point towards carboxyl-end) # TODO test if this works
    if((Symbol(:N_TERMINAL) ∈ filtered_fragment_list[1].flags) 
        && (Symbol(:C_TERMINAL) ∈ filtered_fragment_list[end].flags))
        n_to_c = true
    elseif((Symbol(:N_TERMINAL) ∈ filtered_fragment_list[end].flags) 
        && (Symbol(:C_TERMINAL) ∈ filtered_fragment_list[1].flags))
        n_to_c = false
    else
        throw(ErrorException("c and n terminal are not included in flags (chain $chain_num)"))
    end
    log_info(misc, "n_to_c: $n_to_c")
    
    # for each sheet: find the end and save the index
    prev_idx = nothing
    prev_ss = nothing
    arrow_fragment_indices = []

    for current_fragment in (n_to_c ? fragment_list : reverse(fragment_list))
        current_ss = current_fragment.properties[:SS]
        if(prev_ss!==nothing)
            if(n_to_c 
                && prev_ss==BiochemicalAlgorithms.SecondaryStructure.SHEET 
                && current_ss!=BiochemicalAlgorithms.SecondaryStructure.SHEET)
                push!(arrow_fragment_indices, prev_idx)
            end

            if(!n_to_c 
                && prev_ss!=BiochemicalAlgorithms.SecondaryStructure.SHEET 
                && current_ss==BiochemicalAlgorithms.SecondaryStructure.SHEET)
                push!(arrow_fragment_indices, current_fragment.idx)
            end
        end
        prev_idx = current_fragment.idx
        prev_ss = current_ss
    end

    # find corresponding frames and store the resulting width of the rectangle at that frame
    frames_in_residue_count = 0
    prev_res_idx = sample_to_residue_indices[1]

    a = 1
    while(a<=length(sample_to_residue_indices))
        res_idx = sample_to_residue_indices[a]
        if(res_idx!=prev_res_idx || a==length(sample_to_residue_indices)) # end of the previous residue
            if(a==length(sample_to_residue_indices))
                frames_in_residue_count += 1
            end
            if(prev_res_idx ∈ arrow_fragment_indices)
                num_uniform = Int(round(frames_in_residue_count/3))
                num_arrow = frames_in_residue_count - num_uniform
                uniforms = repeat([1], num_uniform)
                arrow = collect(range(1.5, 0, num_arrow))
                if(n_to_c)
                    append!(arrow_frame_indices, length(rectangle_widths)+num_uniform+1:length(rectangle_widths)+num_uniform+num_arrow)
                    append!(rectangle_widths, uniforms, arrow)
                else
                    append!(arrow_frame_indices, length(rectangle_widths)+1:length(rectangle_widths)+num_arrow)
                    append!(rectangle_widths, reverse!(arrow), uniforms)
                end
                
                # add geometry between uniform and arrow part
                push!(arrow_starts, length(rectangle_widths)-num_arrow+1)
            else
                append!(rectangle_widths, repeat([1], frames_in_residue_count))
            end

            # start new residue
            frames_in_residue_count = 0
            prev_res_idx = res_idx
        end
        frames_in_residue_count += 1
        a += 1
    end

    return rectangle_widths, arrow_starts, arrow_frame_indices
end

function generate_geometry_at_point(
    point::AbstractVector{T}, 
    normal::AbstractVector{T}, 
    binormal::AbstractVector{T}, 
    residue_info_dict::Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}},
    linked_residue_idx::Union{Nothing,Int}, 
    frame_config::Union{Nothing,Tuple{Bool, Bool, Int, Int}}, 
    rectangle_width::T, 
    fixed_color::Union{NTuple{3, Int}, Nothing}, 
    color_dict::Union{Nothing,  Dict{BiochemicalAlgorithms.SecondaryStructure.T, NTuple{3, Int}}, Dict{String, NTuple{3, Int}}},
    config::BackboneConfig) where T
    # generate cross-section vertices
    if(config.backbone_type==BackboneType.BACKBONE)
        circle_points = create_circle_in_local_frame(point, normal, binormal, config.resolution, config.stick_radius)
    elseif(config.backbone_type==BackboneType.RIBBON)
        circle_points = create_ellipse_in_local_frame(point, normal, binormal, config.resolution, T(3)*config.stick_radius, config.stick_radius)
    elseif(config.backbone_type==BackboneType.CARTOON)
        shortcut = false
        if(linked_residue_idx===nothing)
            if((frame_config[1] == frame_config[2]))
                circle_points = stack(repeat([point], config.resolution)) # TODO only 1 instead of resolution vertices would suffice, but then connect_circles_to_tube has to be modified
                shortcut = true
            else
                residue_idx = frame_config[3]
            end
        else
            residue_idx = linked_residue_idx
        end

        if(!shortcut)
            structure = residue_info_dict[residue_idx][2]
            if(structure==BiochemicalAlgorithms.SecondaryStructure.NONE)
                circle_points = create_circle_in_local_frame(point, normal, binormal, config.resolution, config.stick_radius)
            elseif(structure==BiochemicalAlgorithms.SecondaryStructure.HELIX)
                circle_points = create_ellipse_in_local_frame(point, normal, binormal, config.resolution, T(3)*config.stick_radius, T(1.5)*config.stick_radius)
            elseif(structure==BiochemicalAlgorithms.SecondaryStructure.SHEET)
                circle_points = create_rectangle_in_local_frame(point, normal, binormal, config.resolution, T(3)*config.stick_radius * rectangle_width, T(0.5)*config.stick_radius)
            end
        end
    end 
    # add edges
    circle = PlainNonStdMesh(circle_points, Vector{Vector{Int}}(), Vector{NTuple{3, Int}}(undef, config.resolution))

    # color
    if(fixed_color!==nothing)
        color!(circle, fixed_color)
    elseif(config.color==Color.SECONDARY_STRUCTURE)
        if(linked_residue_idx===nothing)
            residue_idx = frame_config[4]
        else
            residue_idx = linked_residue_idx
        end
        structure = residue_info_dict[residue_idx][2]
        color!(circle, color_dict[structure])
    elseif(config.color==Color.RESIDUE)
        if(linked_residue_idx===nothing)
            residue_idx = frame_config[4]
        else
            residue_idx = linked_residue_idx
        end
        color!(circle, color_dict[residue_info_dict[residue_idx][1]])
    end
    return circle
end

# assumes checked config (consistent and adjusted to T)
function prepare_backbone_model(chain::Chain{T}, config::BackboneConfig, fixed_color::Union{Nothing, NTuple{3, Int}} = nothing) where {T<:Real}
    fragment_list = fragments(chain) #TODO fragments/eachfragment

    # construct spline
    local spline
    if(config.spline==Spline.CATMULL_ROM)
        spline = CatmullRom(chain, config.control_point_strategy)
    elseif(config.spline==Spline.CUBIC_B)
        spline = CubicB(chain, config.control_point_strategy)
    end

    vertices_per_unit = T(0.4 * config.resolution / (2*π*config.stick_radius))
    # sample along spline
    spline_points, sample_to_residue_indices::Vector{Union{Int, Nothing}} = calculate_points(spline, vertices_per_unit) #alloc
    velocities = calculate_velocities(spline, vertices_per_unit)
    
    # construct local frames
    local q::Matrix{T}
    local r::Matrix{T}
    local s::Matrix{T}
    if(config.frame==Frame.RMF)
        q, r, s = rmf(spline_points, velocities)
    elseif(config.frame==Frame.SECOND_SPLINE)
        second_spline_points = calculate_minor_points(spline, vertices_per_unit)
        q, r, s = frames_from_two_splines(spline_points, velocities, second_spline_points)


        # sphere_radius = 0.2
        # sphere_mesh = discretize(Sphere{3, Float64}((0,0,0), sphere_radius), RegularDiscretization(6))
        # debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (200, 200, 200)), eachcol(spline_points)))
        # export_mesh_to_ply("main_points_new.ply", debug_mesh)

        # debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (200, 0, 0)), eachcol(second_spline_points)))
        # export_mesh_to_ply("outer_points_new.ply", debug_mesh)
    end
    if(config.color==Color.RAINBOW)
        color_range = range(HSV(0,1,1), stop=HSV(360,1,1), length=size(spline_points, 2))
        rgb_colors = map(hsv-> convert(RGB, hsv), color_range)
        rainbow_colors = map(rgb->map(channel->Int(floor(channel*255)), (rgb.r, rgb.g, rgb.b)), rgb_colors)
    end

    fixed_indices::Vector{Int} = [] # collection of all frames that should not be removed by filtering

    # preparation for arrows
    if(config.backbone_type==BackboneType.CARTOON)
        rectangle_widths, arrow_starts, arrow_frame_indices = compute_frame_widths(fragment_list, sample_to_residue_indices)
        append!(fixed_indices, arrow_frame_indices)
        sort!(fixed_indices)

        for i=eachindex(arrow_starts)
            insertion_idx = arrow_starts[i]
            spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = insert_frame(insertion_idx, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, (config.color==Color.RAINBOW ? rainbow_colors : nothing))
            rectangle_widths[arrow_starts[i]] = T(1.0)

            # adjust indices after the insertion site
            adjust_indices!(fixed_indices, arrow_starts[i])
            adjust_indices!(arrow_starts, arrow_starts[i])

            # newly inserted frame should not be removed by filtering
            insert_sorted!(fixed_indices, insertion_idx)

        end

    end

    # add frames when secondary structure changes # TODO what if n_to_c is false?
    if(config.backbone_type==BackboneType.CARTOON)
        frame_config = Dict{Int, Tuple{Bool, Bool, Int, Int}}()
        prev_res_idx = sample_to_residue_indices[1]
        a = 1
        while(a<=length(sample_to_residue_indices))
            res_idx = sample_to_residue_indices[a]
            prev_ss = spline.residue_info_dict[prev_res_idx][2]
            curr_ss = spline.residue_info_dict[res_idx][2]

            if(res_idx!=prev_res_idx && prev_ss!=curr_ss)
                ss_a = prev_ss
                ss_b = curr_ss
                small_to_large = ss_a==BiochemicalAlgorithms.SecondaryStructure.NONE || ss_a==BiochemicalAlgorithms.SecondaryStructure.SHEET
                #log_info(extra_frames, "$ss_a -> $ss_b, small_to_large: $small_to_large")

                if(small_to_large)
                    insertion_idx = a-1
                    #log_info(extra_frames, insertion_idx, fixed_indices)
                    spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = insert_frame(insertion_idx, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, (config.color==Color.RAINBOW ? rainbow_colors : nothing))
                    spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = insert_frame(insertion_idx, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, (config.color==Color.RAINBOW ? rainbow_colors : nothing))
                    sample_to_residue_indices[insertion_idx+1] = nothing
                    sample_to_residue_indices[insertion_idx+2] = nothing
                    
                    frame_config[insertion_idx+1] = (small_to_large, true, prev_res_idx, res_idx)
                    frame_config[insertion_idx+2] = (small_to_large, false, res_idx, res_idx)

                    adjust_indices!(fixed_indices, insertion_idx+1)
                    adjust_indices!(fixed_indices, insertion_idx+1)

                    insert_sorted!(fixed_indices, insertion_idx)
                    insert_sorted!(fixed_indices, insertion_idx+1)
                    insert_sorted!(fixed_indices, insertion_idx+2)
                    #log_info(extra_frames, insertion_idx, fixed_indices)
                else
                    insertion_idx = a
                    #log_info(extra_frames, insertion_idx, fixed_indices)
                    spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = insert_frame(insertion_idx, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, (config.color==Color.RAINBOW ? rainbow_colors : nothing))
                    spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = insert_frame(insertion_idx, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, (config.color==Color.RAINBOW ? rainbow_colors : nothing))
                    sample_to_residue_indices[insertion_idx] = nothing
                    sample_to_residue_indices[insertion_idx+1] = nothing
                    
                    frame_config[insertion_idx] = (small_to_large, true, prev_res_idx, prev_res_idx)
                    frame_config[insertion_idx+1] = (small_to_large, false, res_idx, prev_res_idx)

                    adjust_indices!(fixed_indices, insertion_idx)
                    adjust_indices!(fixed_indices, insertion_idx)

                    insert_sorted!(fixed_indices, insertion_idx)
                    insert_sorted!(fixed_indices, insertion_idx+1)
                    insert_sorted!(fixed_indices, insertion_idx+2)
                    #log_info(extra_frames, insertion_idx, fixed_indices)
                    
                end
                a+=2
                #log_info(extra_frames)
            
            end
            prev_res_idx = res_idx
            a+=1
        end


    end


    # filter
    if(config.filter==Filter.ANGLE)
        # create a list of points that are fixed and cannot be removed by the filter
        insert!(fixed_indices, 1, 1) # begin and end cannot be dropped
        push!(fixed_indices, size(spline_points, 2))

        if(config.color==Color.SECONDARY_STRUCTURE && config.backbone_type!=BackboneType.CARTOON)
            prev_res_idx = sample_to_residue_indices[1]
            prev_ss = spline.residue_info_dict[sample_to_residue_indices[1]][2]
            for (i, res_idx) in enumerate(sample_to_residue_indices)

                if(res_idx!=prev_res_idx && spline.residue_info_dict[res_idx][2]!=prev_ss)
                    push!(fixed_indices, i-1, i) # Color changes at secondary structure changes. 
                    prev_ss = spline.residue_info_dict[res_idx][2]
                    prev_res_idx = res_idx
                end
            end
        end
        if(config.color==Color.RESIDUE)
            prev_res_idx = sample_to_residue_indices[1]
            prev_ss = spline.residue_info_dict[sample_to_residue_indices[1]][2]

            prev_was_nothing = false
            for (i, resIdx) in enumerate(sample_to_residue_indices)

                if(resIdx===nothing)
                    prev_was_nothing = true
                    continue
                end


                if(resIdx!=prev_res_idx && !prev_was_nothing)
                    push!(fixed_indices, i-1, i) # Color changes at residue boundaries. The color should not be interpolated over a larger distance.
                    prev_res_idx = resIdx
                end

                prev_was_nothing = false
            end
        end

        local remaining_indices::Set{Int}
        if(config.color==Color.RAINBOW)
            remaining_indices = filter_points_threshold(q, r, fixed_indices, rainbow_colors) # prevents too large distances in colors as well
        else
            remaining_indices = filter_points_threshold(q, r, fixed_indices)
        end
        # log_info(point_filter, "Remaining points: $(length(remaining_indices))/$(size(spline_points, 2))\n Filtered: $(size(spline_points, 2)-length(remaining_indices)) ($(length(fixed_indices)) fixed)")
        # spline_points = spline_points[:, remaining_indices]
        # sample_to_residue_indices = sample_to_residue_indices[remaining_indices]
        # q = q[:, remaining_indices]
        # r = r[:, remaining_indices]
        # s = s[:, remaining_indices]
        # rectangle_widths = rectangle_widths[remaining_indices]
        # if(config.color==Color.RAINBOW)
        #     rainbow_colors = rainbow_colors[remaining_indices]
        # end
    end

    color_dict = nothing
    if(config.color==Color.SECONDARY_STRUCTURE)
        color_dict = get_structure_color_mapping()
    elseif(config.color==Color.RESIDUE)
        color_dict = get_amino_acid_color_mapping()
    end


    log_info(types, "Type of spline points: ", typeof(spline_points))

    circles::Vector{PlainNonStdMesh{T}} = []

    # framesA::Vector{ColoredMesh} = []
    # framesB::Vector{ColoredMesh} = []
    # push!(framesB, local_frame_mesh(spline_points[:, 1], q[:, 1], r[:, 1], s[:, 1]))

    # iterate and create vertices
    for current_index=axes(spline_points, 2)
        if(config.filter!=Filter.NONE && current_index ∉ remaining_indices)
            continue
        end
        # sanity check: frame should be orthogonal
        @views if(!approx_zero(dot(q[:, current_index], r[:, current_index])) || !approx_zero(dot(q[:, current_index], s[:, current_index])) || !approx_zero(dot(s[:, current_index], r[:, current_index])))
            log_info(frame_rotation, current_index, " wrong angles ", dot(q[:, current_index], r[:, current_index]), " ", dot(q[:, current_index], s[:, current_index]), " ", dot(s[:, current_index], r[:, current_index]), " # ", q[:, current_index], " ", r[:, current_index], " ", s[:, current_index])
        end

        # debug output
        # frameA = local_frame_mesh(spline_points[:, current_index], q[:, current_index], r[:, current_index], s[:, current_index])
        # push!(framesA, frameA)

        
        if(config.color==Color.RAINBOW)
            fixed_color = rainbow_colors[current_index]
        elseif(fixed_color===nothing && (config.color==Color.UNIFORM || config.color==Color.CHAIN))
            fixed_color = (0, 0, 255)
        end

        circle = @views generate_geometry_at_point(
            spline_points[:, current_index], 
            r[:, current_index], 
            s[:, current_index], 
            spline.residue_info_dict,
            sample_to_residue_indices[current_index],
            (config.backbone_type==BackboneType.CARTOON && current_index ∈ keys(frame_config)) ? frame_config[current_index] : nothing,
            config.backbone_type==BackboneType.CARTOON ? rectangle_widths[current_index] : T(1.0),
            fixed_color, 
            color_dict,
            config)

        push!(circles, circle)
    end

    spline_mesh = connect_circles_to_tube(circles, @views (spline_points[:, 1], spline_points[:, end]))
    log_info(types, "Type of spline mesh: ", typeof(spline_mesh))

    # ----- debug export -----
    # cs = merge_multiple_meshes(circle_meshes)
    # export_mesh_representation_to_ply("circles.ply", Representation(cs))

    # fs = reduce(merge, framesA)
    # export_mesh_to_ply("framesA.ply", fs)

    # fs = reduce(merge, framesB)
    # export_mesh_to_ply("framesB.ply", fs)
    # ------------------------
    
    return spline_mesh
end

function prepare_backbone_model(
    ac::System{T}, config::BackboneConfig) where {T<:Real}

    start_time = now()
    config.stick_radius = T(config.stick_radius)


    log_info(types, "Type: ", T)

    if(config.color==Color.ELEMENT)
        throw(ArgumentError("backbone-based models cannot be colored by elements of individual atoms"))
    end
    if(config.control_point_strategy==ControlPoints.C_ALPHA && config.frame==Frame.SECOND_SPLINE)
        throw(ArgumentError("for a second spline, ControlPoints.MID_POINTS is mandatory"))
    end

    

    if(config.color==Color.UNIFORM)
        uniform_color = (255, 0, 0)
    elseif(config.color==Color.CHAIN)
        chain_colors = map(c->map(channel->Int(channel*255), (c.r, c.g, c.b)), collect(distinguishable_colors(nchains(ac)+1))[2:end]) #alloc
    end

    chain_meshes::Vector{PlainMesh{T}} = []
    for (chain_num, chain) in enumerate(BiochemicalAlgorithms.chains(ac))
        try
            color = nothing
            if(config.color==Color.UNIFORM)
                color = uniform_color
            elseif(config.color==Color.CHAIN)
                color = chain_colors[chain_num]
            end

            chain_mesh = prepare_backbone_model(chain, config, color)

            push!(chain_meshes, chain_mesh)
        catch e
            if e isa ErrorException
                log_warning("Skipped chain $(chain.name), because an error occured: $(e.msg)")
            else
                rethrow(e)
            end
        end
    end
    if(length(chain_meshes)==0)
        throw(ErrorException("No chain meshes were generated. "))
    end
    temp = merge_multiple_meshes(chain_meshes)
    log_info(types, "Type of result: ", typeof(temp))

    log_info(time_info, "Generated backbone mesh in $((now()-start_time).value/1000) seconds. ($(size(temp.vertices, 2)) vertices)")

    return temp

end