export prepare_backbone_model

default_config = BackboneConfig(0.2, 
12, 
BackboneType.CARTOON, 
Color.SECONDARY_STRUCTURE, 
Spline.CUBIC_B, 
ControlPoints.MID_POINTS, 
Frame.SECOND_SPLINE, 
Filter.ANGLE)

function insert_sorted!(array, elem)
    index = searchsortedfirst(array, elem)
    insert!(array, index, elem)
end

function insert_frame(index, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = nothing)
    if(index<1 || index>size(spline_points, 2))
        throw(ArgumentError("index has to be in [1, $(size(spline_points, 2))], but was $index"))
    end
    spline_points = hcat(spline_points[:, 1:index], spline_points[:, index], spline_points[:, index+1:end])
    q = hcat(q[:, 1:index], q[:, index], q[:, index+1:end])
    r = hcat(r[:, 1:index], r[:, index], r[:, index+1:end])
    s = hcat(s[:, 1:index], s[:, index], s[:, index+1:end])
    rectangle_widths = vcat(rectangle_widths[1:index], rectangle_widths[index], rectangle_widths[index+1:end])
    sample_to_residue_indices = vcat(sample_to_residue_indices[1:index], sample_to_residue_indices[index], sample_to_residue_indices[index+1:end])
    if(rainbow_colors !== nothing)
        rainbow_colors = vcat(rainbow_colors[1:index], rainbow_colors[index], rainbow_colors[index+1:end])
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

function compute_frame_widths(fragment_list, sample_to_residue_indices)
    rectangle_widths = []
    arrow_starts = []
    arrow_frame_indices = []

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
    fragment_indices = collect(1:length(fragment_list))
    if(!n_to_c)
        reverse!(fragment_indices)
    end
    prev_i = nothing
    prev_ss = nothing
    arrow_fragment_indices = []
    for i in fragment_indices
        current_ss = fragment_list[i].properties[:SS]
        if(prev_ss!==nothing)
            if(n_to_c 
                && prev_ss==BiochemicalAlgorithms.SecondaryStructure.SHEET 
                && current_ss!=BiochemicalAlgorithms.SecondaryStructure.SHEET)
                push!(arrow_fragment_indices, prev_i)
            end

            if(!n_to_c 
                && prev_ss!=BiochemicalAlgorithms.SecondaryStructure.SHEET 
                && current_ss==BiochemicalAlgorithms.SecondaryStructure.SHEET)
                push!(arrow_fragment_indices, i)
            end
        end
        prev_i = i
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

function prepare_backbone_model(
    ac::AbstractAtomContainer{T}, config::BackboneConfig=default_config) where {T<:Real}

    start_time = now()
    U = Float64
    if(T <: AbstractFloat)
        U = T
    end

    log_info(types, "Types: ", T, " ", U)

    if(config.color==Color.ELEMENT)
        throw(ArgumentError("backbone-based models cannot be colored by elements of individual atoms"))
    end
    if(config.control_point_strategy==ControlPoints.C_ALPHA && config.frame==Frame.SECOND_SPLINE)
        throw(ArgumentError("for a second spline, ControlPoints.MID_POINTS is mandatory"))
    end

    vertices_per_unit = 0.4 * config.resolution / (2*π*config.stick_radius)
    # circle_prototype = discretize(Sphere(Point(U(0),U(0)), stick_radius), RegularDiscretization(resolution))
    # circle_mesh = PlainNonStdMesh(lift_into_3d(circle_prototype))
    # cap_mesh = PlainMesh(Hemisphere(stick_radius, resolution, U))

    #log_info(types, "Types of proto meshes: ", typeof(circle_mesh), " ", typeof(cap_mesh))
    if(config.color==Color.UNIFORM)
        uniform_color = (255, 0, 0)
    elseif(config.color==Color.CHAIN)
        chain_colors = map(c->map(channel->Int(channel*255), (c.r, c.g, c.b)), collect(distinguishable_colors(nchains(ac)+1))[2:end])
    elseif(config.color==Color.SECONDARY_STRUCTURE)
        structure_color_mapping = get_structure_color_mapping()
    elseif(config.color==Color.RESIDUE)
        amino_acid_color_mapping = get_amino_acid_color_mapping()
    end

    chain_meshes::Vector{PlainMesh{U}} = []
    aborted_chain_count = 0
    for (chain_num, chain) in enumerate(BiochemicalAlgorithms.chains(ac))
        fragment_list = fragments(chain)

        # construct spline
        local spline
        try
            if(config.spline==Spline.CATMULL_ROM)
                spline = CatmullRom(chain, config.control_point_strategy)
            elseif(config.spline==Spline.CUBIC_B)
                spline = CubicB(chain, config.control_point_strategy)
            end
        catch e
            aborted_chain_count += 1
            continue
        end


        # sample along spline
        spline_points::Matrix{U}, sample_to_residue_indices::Vector{Union{Int, Nothing}} = calculate_points(spline, vertices_per_unit)
        velocities::Matrix{U} = calculate_velocities(spline, vertices_per_unit)
        
        # construct local frames
        local q::Matrix{U}
        local r::Matrix{U}
        local s::Matrix{U}
        if(config.frame==Frame.RMF)
            q, r, s = rmf(spline_points, velocities)
        elseif(config.frame==Frame.SECOND_SPLINE)
            second_spline_points::Matrix{U} = calculate_minor_points(spline, vertices_per_unit)
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

        fixed_indices = [] # collection of all frames that should not be removed by filtering

        # preparation for arrows
        if(config.backbone_type==BackboneType.CARTOON)
            rectangle_widths, arrow_starts, arrow_frame_indices = compute_frame_widths(fragment_list, sample_to_residue_indices)
            append!(fixed_indices, arrow_frame_indices)
            sort!(fixed_indices)

            for i=eachindex(arrow_starts)
                insertion_idx = arrow_starts[i]
                spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = insert_frame(insertion_idx, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, (config.color==Color.RAINBOW ? rainbow_colors : nothing))
                rectangle_widths[arrow_starts[i]] = 1.0

                # adjust indices after the insertion site
                adjust_indices!(fixed_indices, arrow_starts[i])
                adjust_indices!(arrow_starts, arrow_starts[i])

                # newly inserted frame should not be removed by filtering
                insert_sorted!(fixed_indices, insertion_idx)

            end

        end

        # add frames when secondary structure changes # TODO what if n_to_c is false?
        if(config.backbone_type==BackboneType.CARTOON)
            frame_config = Dict()
            prev_res_idx = sample_to_residue_indices[1]
            a = 1
            while(a<=length(sample_to_residue_indices))
                res_idx = sample_to_residue_indices[a]
                if(res_idx!=prev_res_idx && fragment_list[prev_res_idx].properties[:SS]!=fragment_list[res_idx].properties[:SS])
                    ss_a = fragment_list[prev_res_idx].properties[:SS]
                    ss_b = fragment_list[res_idx].properties[:SS]
                    small_to_large = ss_a==BiochemicalAlgorithms.SecondaryStructure.NONE || ss_a==BiochemicalAlgorithms.SecondaryStructure.SHEET
                    #log_info(extra_frames, "$ss_a -> $ss_b, small_to_large: $small_to_large")

                    if(small_to_large)
                        insertion_idx = a-1
                        #log_info(extra_frames, insertion_idx, fixed_indices)
                        spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = insert_frame(insertion_idx, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, (config.color==Color.RAINBOW ? rainbow_colors : nothing))
                        spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = insert_frame(insertion_idx, spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, (config.color==Color.RAINBOW ? rainbow_colors : nothing))
                        sample_to_residue_indices[insertion_idx+1] = nothing
                        sample_to_residue_indices[insertion_idx+2] = nothing
                        
                        frame_config[insertion_idx+1] = (small_to_large, true, fragment_list[prev_res_idx], fragment_list[res_idx])
                        frame_config[insertion_idx+2] = (small_to_large, false, fragment_list[res_idx], fragment_list[res_idx])

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
                        
                        frame_config[insertion_idx] = (small_to_large, true, fragment_list[prev_res_idx], fragment_list[prev_res_idx])
                        frame_config[insertion_idx+1] = (small_to_large, false, fragment_list[res_idx], fragment_list[prev_res_idx])

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
                prev_ss = fragment_list[sample_to_residue_indices[1]].properties[:SS]
                for (i, res_idx) in enumerate(sample_to_residue_indices)

                    if(res_idx!=prev_res_idx && fragment_list[res_idx].properties[:SS]!=prev_ss)
                        push!(fixed_indices, i-1, i) # Color changes at secondary structure changes. 
                        prev_ss = fragment_list[res_idx].properties[:SS]
                        prev_res_idx = res_idx
                    end
                end
            end
            if(config.color==Color.RESIDUE)
                prev_res_idx = sample_to_residue_indices[1]
                prev_ss = fragment_list[sample_to_residue_indices[1]].properties[:SS]

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
            sort!(fixed_indices)
            unique!(fixed_indices)

            local remaining_indices::Vector{Int}
            if(config.color==Color.RAINBOW)
                remaining_indices = filter_points_threshold(spline_points, q, r, s, fixed_indices, rainbow_colors) # prevents too large distances in colors as well
            else
                remaining_indices = filter_points_threshold(spline_points, q, r, s, fixed_indices)
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


        log_info(types, "Type of spline points: ", typeof(spline_points))

        circles::Vector{PlainNonStdMesh{U}} = []

        # framesA::Vector{ColoredMesh} = []
        # framesB::Vector{ColoredMesh} = []
        # push!(framesB, local_frame_mesh(spline_points[:, 1], q[:, 1], r[:, 1], s[:, 1]))

        # iterate and create vertices
        for current_index=axes(spline_points, 2)
            if(config.filter!=Filter.NONE && current_index ∉ remaining_indices)
                continue
            end
            # sanity check: frame should be orthogonal
            if(!approx_zero(dot(q[:, current_index], r[:, current_index])) || !approx_zero(dot(q[:, current_index], s[:, current_index])) || !approx_zero(dot(s[:, current_index], r[:, current_index])))
                log_info(frame_rotation, current_index, " wrong angles ", dot(q[:, current_index], r[:, current_index]), " ", dot(q[:, current_index], s[:, current_index]), " ", dot(s[:, current_index], r[:, current_index]), " # ", q[:, current_index], " ", r[:, current_index], " ", s[:, current_index])
            end

            # debug output
            # frameA = local_frame_mesh(spline_points[:, current_index], q[:, current_index], r[:, current_index], s[:, current_index])
            # push!(framesA, frameA)

            # generate cross-section
            if(config.backbone_type==BackboneType.BACKBONE)
                circle_points = create_circle_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, config.stick_radius)
            elseif(config.backbone_type==BackboneType.RIBBON)
                circle_points = create_ellipse_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, T(3)*config.stick_radius, config.stick_radius)
            elseif(config.backbone_type==BackboneType.CARTOON)
                shortcut = false
                if(sample_to_residue_indices[current_index]===nothing)
                    if((frame_config[current_index][1] == frame_config[current_index][2]))
                        circle_points = stack(repeat([spline_points[:, current_index]], config.resolution)) # TODO only 1 instead of resolution vertices would suffice, but then connect_circles_to_tube has to be modified
                        shortcut = true
                    else
                        residue = frame_config[current_index][3]
                    end
                else
                    residue = fragment_list[sample_to_residue_indices[current_index]]
                end

                if(!shortcut)
                    structure = residue.properties[:SS]
                    if(structure==BiochemicalAlgorithms.SecondaryStructure.NONE)
                        circle_points = create_circle_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, config.stick_radius)
                    elseif(structure==BiochemicalAlgorithms.SecondaryStructure.HELIX)
                        circle_points = create_ellipse_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, T(3)*config.stick_radius, T(1.5)*config.stick_radius)
                    elseif(structure==BiochemicalAlgorithms.SecondaryStructure.SHEET)
                        circle_points = create_rectangle_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, T(3)*config.stick_radius * T(rectangle_widths[current_index]), T(0.5)*config.stick_radius)
                    end
                end
            end 
            circle = PlainNonStdMesh(circle_points, Vector{Vector{Int}}(), Vector{NTuple{3, Int}}())

            if(config.color==Color.UNIFORM)
                color!(circle, uniform_color)
            elseif(config.color==Color.CHAIN)
                color!(circle, chain_colors[chain_num])
            elseif(config.color==Color.RAINBOW)
                color!(circle, rainbow_colors[current_index])
            elseif(config.color==Color.SECONDARY_STRUCTURE)
                if(sample_to_residue_indices[current_index]===nothing)
                    residue = frame_config[current_index][4]
                else
                    residue = fragment_list[sample_to_residue_indices[current_index]]
                end
                structure = residue.properties[:SS]
                color!(circle, structure_color_mapping[structure])
            elseif(config.color==Color.RESIDUE)
                if(sample_to_residue_indices[current_index]===nothing)
                    residue = frame_config[current_index][4]
                else
                    residue = fragment_list[sample_to_residue_indices[current_index]]
                end
                color!(circle, amino_acid_color_mapping[residue.name])
            end

            push!(circles, circle)
        end

        spline_mesh = connect_circles_to_tube(circles, [spline_points[:, 1], spline_points[:, end]])
        log_info(types, "Type of spline mesh: ", typeof(spline_mesh))

        # ----- debug export -----
        # cs = merge_multiple_meshes(circle_meshes)
        # export_mesh_representation_to_ply("circles.ply", Representation(cs))

        # fs = reduce(merge, framesA)
        # export_mesh_to_ply("framesA.ply", fs)

        # fs = reduce(merge, framesB)
        # export_mesh_to_ply("framesB.ply", fs)
        # ------------------------
        
        push!(chain_meshes, spline_mesh)
    end
    if(length(chain_meshes)==0)
        throw(ErrorException("all chains had too few c_alpha atoms. Mesh could not be generated. "))
    end
    if(aborted_chain_count>0)
        log_warning("Skipped $aborted_chain_count chains because they did not conain enough c_alpha atoms. ")
    end
    temp = merge_multiple_meshes(chain_meshes)
    result = Representation(temp)
    log_info(types, "Type of result: ", typeof(result))

    log_info(time_info, "Generated backbone mesh in $((now()-start_time).value/1000) seconds. ($(length(result.vertices) ÷ 3) vertices)")

    return result

end