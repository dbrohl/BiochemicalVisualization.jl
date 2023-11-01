export prepare_backbone_model

default_config = BackboneConfig(0.2, 
12, 
BackboneType.CARTOON, 
Color.SECONDARY_STRUCTURE, 
Spline.CUBIC_B, 
ControlPoints.MID_POINTS, 
Frame.SECOND_SPLINE, 
Filter.ANGLE)

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
                    append!(arrow_frame_indices, length(rectangle_widths)+num_uniform:length(rectangle_widths)+num_uniform+num_arrow)
                    append!(rectangle_widths, uniforms, arrow)
                else
                    append!(arrow_frame_indices, length(rectangle_widths)+1:length(rectangle_widths)+num_arrow+1)
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
    for (chain_num, chain) in enumerate(BiochemicalAlgorithms.chains(ac))
        println()
        fragment_list = fragments(chain)
        # TODO zu kurze chains abfangen

        # construct spline
        if(config.spline==Spline.CATMULL_ROM)
            spline = CatmullRom(chain, config.control_point_strategy)
        elseif(config.spline==Spline.CUBIC_B)
            spline = CubicB(chain, config.control_point_strategy)
        end

        # sample along spline
        spline_points::Matrix{U}, sample_to_residue_indices::Vector{Int} = calculate_points(spline, vertices_per_unit)
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
                spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, rainbow_colors = insert_frame(arrow_starts[i], spline_points, q, r, s, rectangle_widths, sample_to_residue_indices, (config.color==Color.RAINBOW ? rainbow_colors : nothing))
                rectangle_widths[arrow_starts[i]] = 1.0

                # adjust indices after the insertion site
                a = 1
                inserted = false
                while a<length(fixed_indices)
                    if(fixed_indices[a]>=arrow_starts[i])
                        if(!inserted)
                            insert!(fixed_indices, a, arrow_starts[i])
                            inserted = true
                            a += 1
                        end
                        fixed_indices[a] += 1
                    end
                    a += 1
                end

                for j=eachindex(arrow_starts)
                    if(arrow_starts[j]> arrow_starts[i])
                        arrow_starts[j] += 1
                    end
                end
            end
        end

        # TODO add frames when secondary structure changes


        # filter
        if(config.filter==Filter.ANGLE)
            # create a list of points that are fixed and cannot be removed by the filter # TODO in 2lwq, one sheet has a too short residue coloring
            push!(fixed_indices, 1, size(spline_points, 2)) # begin and end cannot be dropped
            if(config.backbone_type==BackboneType.CARTOON || config.color==Color.SECONDARY_STRUCTURE || config.color==Color.RESIDUE)
                prev_res_idx = sample_to_residue_indices[1]
                prev_ss = fragment_list[sample_to_residue_indices[1]].properties[:SS]
                for (i, resIdx) in enumerate(sample_to_residue_indices)


                    if(resIdx!=prev_res_idx)
                        residue = fragment_list[resIdx]
                        if(config.color==Color.RESIDUE)
                            push!(fixed_indices, i-1, i) # Color changes at residue boundaries. The color should not be interpolated over a larger distance.
                        end
                        if((config.backbone_type==BackboneType.CARTOON || config.color==Color.SECONDARY_STRUCTURE)
                            && residue.properties[:SS]!=prev_ss)
                            push!(fixed_indices, i-1, i) # Color or cross section changes at secondary structure changes. 
                            prev_ss = residue.properties[:SS]
                        end
                        prev_res_idx = resIdx
                    end
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
            log_info(point_filter, "Remaining points: $(length(remaining_indices))/$(size(spline_points, 2))\n Filtered: $(size(spline_points, 2)-length(remaining_indices)) ($(length(fixed_indices)) fixed)")
            spline_points = spline_points[:, remaining_indices]
            sample_to_residue_indices = sample_to_residue_indices[remaining_indices]
            q = q[:, remaining_indices]
            r = r[:, remaining_indices]
            s = s[:, remaining_indices]
            rectangle_widths = rectangle_widths[remaining_indices]
            if(config.color==Color.RAINBOW)
                rainbow_colors = rainbow_colors[remaining_indices]
            end
        end


        log_info(types, "Type of spline points: ", typeof(spline_points))

        circles::Vector{PlainNonStdMesh{U}} = []

        # framesA::Vector{ColoredMesh} = []
        # framesB::Vector{ColoredMesh} = []
        # push!(framesB, local_frame_mesh(spline_points[:, 1], q[:, 1], r[:, 1], s[:, 1]))

        # iterate and create vertices
        for current_index=axes(spline_points, 2) # start- and endcap bases are the first and last splinepoints #TODO offset wieder auf 2, wenn Endkappen wieder funktionieren

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
                residue = fragment_list[sample_to_residue_indices[current_index]]
                structure = residue.properties[:SS]
                if(structure==BiochemicalAlgorithms.SecondaryStructure.NONE)
                    circle_points = create_circle_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, config.stick_radius)
                elseif(structure==BiochemicalAlgorithms.SecondaryStructure.HELIX)
                    circle_points = create_ellipse_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, T(3)*config.stick_radius, T(1.5)*config.stick_radius)
                elseif(structure==BiochemicalAlgorithms.SecondaryStructure.SHEET)
                    circle_points = create_rectangle_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, T(3)*config.stick_radius * T(rectangle_widths[current_index]), T(0.5)*config.stick_radius)
                end
            end
            circle = PlainNonStdMesh(circle_points, Vector{Vector{Int}}(), Vector{NTuple{3, Int}}())
            #TODO extra geometry at SS and arrow boundaries?

            if(config.color==Color.UNIFORM)
                color!(circle, uniform_color)
            elseif(config.color==Color.CHAIN)
                color!(circle, chain_colors[chain_num])
            elseif(config.color==Color.RAINBOW)
                color!(circle, rainbow_colors[current_index])
            elseif(config.color==Color.SECONDARY_STRUCTURE)
                residue = fragment_list[sample_to_residue_indices[current_index]]
                color!(circle, structure_color_mapping[residue.properties[:SS]])
            elseif(config.color==Color.RESIDUE)
                residue = fragment_list[sample_to_residue_indices[current_index]]
                color!(circle, amino_acid_color_mapping[residue.name])
            end

            push!(circles, circle)
        end

        spline_mesh = connect_circles_to_tube(circles)
        log_info(types, "Type of spline mesh: ", typeof(spline_mesh))

        # ----- debug export -----
        # cs = merge_multiple_meshes(circle_meshes)
        # export_mesh_representation_to_ply("circles.ply", Representation(cs))

        # fs = reduce(merge, framesA)
        # export_mesh_to_ply("framesA.ply", fs)

        # fs = reduce(merge, framesB)
        # export_mesh_to_ply("framesB.ply", fs)

        # ----- add hemispheres to both ends -----

    #     # position
    #     TODO start_cap = deepcopy(cap_mesh)
    #     rotate_in_direction!(start_cap, -(spline_points[:, 2]-spline_points[:, 1]))
    #     translate!(start_cap, spline_points[:, 1])
    #     end_cap = deepcopy(cap_mesh)
    #     rotate_in_direction!(end_cap, -(spline_points[:, end-1]-spline_points[:, end]))
    #     translate!(end_cap, spline_points[:, end-1])

    #     color!(start_cap, chain_colors[chain_num])
    #     color!(end_cap, chain_colors[chain_num])
    #     # color!(start_cap, (0, 255, 0))
    #     # color!(end_cap, (0, 0, 255))

    #     log_info(types, "Typeof shifted cap: ", typeof(start_cap), " ", typeof(end_cap))

    #     # merge meshes

    #     points = [spline_mesh.vertices start_cap.vertices end_cap.vertices]
    #     colors = [spline_mesh.colors; start_cap.colors; end_cap.colors]

    #     connections::Matrix{Int} = Matrix(undef, 3, nconnections(spline_mesh) + nconnections(start_cap) + nconnections(end_cap) + 4 * resolution)
        
    #     position = 1
    #     offset = 0
    #     connections[:, position : position + nconnections(spline_mesh)-1] = spline_mesh.connections
    #     position += nconnections(spline_mesh)
    #     offset += nvertices(spline_mesh)

    #     connections[:, position : position + nconnections(start_cap)-1] = (start_cap.connections .+ offset)
    #     position += nconnections(start_cap)
    #     offset += nvertices(start_cap)

    #     connections[:, position : position + nconnections(end_cap)-1] = (end_cap.connections .+ offset)
    #     position += nconnections(end_cap)
    #     offset += nvertices(end_cap)

    #     # add new connections

    #     start_cap_base_indices = nvertices(spline_mesh)+1 : nvertices(spline_mesh)+resolution
    #     end_cap_base_indices = nvertices(spline_mesh)+nvertices(start_cap)+1 : nvertices(spline_mesh)+nvertices(start_cap)+resolution
    #     first_circle_indices = 1 : resolution
    #     last_circle_indices = nvertices(spline_mesh)-resolution+1 : nvertices(spline_mesh)

    #     shift, flip = determine_offset(points[:, start_cap_base_indices[1]], points[:, start_cap_base_indices[2]], points[:, first_circle_indices])
    #     first_circle_indices = circshift(first_circle_indices, -shift)
    #     if(flip)
    #         reverse!(first_circle_indices)
    #     end

    #     shift, flip = determine_offset(points[:, end_cap_base_indices[1]], points[:, end_cap_base_indices[2]], points[:, last_circle_indices])
    #     end_cap_base_indices = circshift(end_cap_base_indices, -shift)
    #     if(flip)
    #         reverse!(end_cap_base_indices)
    #     end


    #     connections[1, position:position+resolution-1] = start_cap_base_indices'
    #     connections[2, position:position+resolution-1] = first_circle_indices'
    #     connections[3, position:position+resolution-1] = circshift(first_circle_indices, 1)'
    #     position += resolution

    #     connections[1, position:position+resolution-1] = start_cap_base_indices'
    #     connections[2, position:position+resolution-1] = circshift(first_circle_indices, 1)'
    #     connections[3, position:position+resolution-1] = circshift(start_cap_base_indices, 1)'
    #     position += resolution

    #     connections[1, position:position+resolution-1] = end_cap_base_indices'
    #     connections[2, position:position+resolution-1] = last_circle_indices'
    #     connections[3, position:position+resolution-1] = circshift(last_circle_indices, 1)'
    #     position += resolution

    #     connections[1, position:position+resolution-1] = end_cap_base_indices'
    #     connections[2, position:position+resolution-1] = circshift(last_circle_indices, 1)'
    #     connections[3, position:position+resolution-1] = circshift(end_cap_base_indices, 1)'
    #     position += resolution
        
    #     push!(chain_meshes, PlainMesh(points, connections, colors))
    #     log_info(types, "Type of capped spline mesh: ", typeof(chain_meshes[end]))
        push!(chain_meshes, spline_mesh)
    end
    temp = merge_multiple_meshes(chain_meshes)
    result = Representation(temp)
    log_info(types, "Type of result: ", typeof(result))

    log_info(time_info, "Generated backbone mesh in $((now()-start_time).value/1000) seconds. ($(length(result.vertices) ÷ 3) vertices)")

    return result

end