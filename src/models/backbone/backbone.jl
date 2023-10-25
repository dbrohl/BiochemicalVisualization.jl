export rotation_test, prepare_backbone_model


FLIP_COLOR::NTuple{3, Int} = (200, 200, 200)
BREAK_COLOR::NTuple{3, Int} = (255, 0, 0)
COUNT_COLOR_START::NTuple{3, Int} = (255, 0, 255)
COUNT_COLOR_0::NTuple{3, Int} = (0, 0, 255)
COUNT_COLOR_1::NTuple{3, Int} = (0, 255, 255)
BACKGROUND_COLOR::NTuple{3, Int} = (60, 60, 60)

structure_color_mapping = Dict(
    BiochemicalAlgorithms.SecondaryStructure.NONE => (255, 255, 255), 
    BiochemicalAlgorithms.SecondaryStructure.HELIX => (255, 75, 120), 
    BiochemicalAlgorithms.SecondaryStructure.SHEET => (255, 150, 0))

amino_acid_color_mapping = Dict()
for ((aa, aa_three_letters, aa_one_letter), color) in zip(values(BiochemicalAlgorithms.AminoAcidProperties), distinguishable_colors(length(AminoAcidProperties)))
    amino_acid_color_mapping[aa_three_letters] = map(channel->Int(channel*255), (color.r, color.g, color.b))
end


function circlesToSimpleMesh(circles)
    verts = hcat(map(c->c.vertices, circles)...)
    verts = [Point3(p...) for p in eachcol(verts)]

    i=0
    connects::Vector{Connectivity} = []
    for c in circles
        for connection in c.connections
            push!(connects, connect(Tuple(connection .+ i)))
        end
        i+=size(c.vertices, 2)
    end

    return SimpleMesh(verts, connects)

end



function prepare_backbone_model(
    ac::AbstractAtomContainer{T}, config::BackboneConfig) where {T<:Real}

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

    vertices_per_unit = config.resolution / (2*π*config.stick_radius)
    # circle_prototype = discretize(Sphere(Point(U(0),U(0)), stick_radius), RegularDiscretization(resolution))
    # circle_mesh = PlainNonStdMesh(lift_into_3d(circle_prototype))
    # cap_mesh = PlainMesh(Hemisphere(stick_radius, resolution, U))

    #log_info(types, "Types of proto meshes: ", typeof(circle_mesh), " ", typeof(cap_mesh))
    if(config.color==Color.UNIFORM)
        uniform_color = (255, 0, 0)
    elseif(config.color==Color.CHAIN)
        chain_colors = map(c->map(channel->Int(channel*255), (c.r, c.g, c.b)), collect(distinguishable_colors(nchains(ac)+1))[2:end])
    end

    chain_meshes::Vector{PlainMesh{U}} = []
    for (chain_num, chain) in enumerate(BiochemicalAlgorithms.chains(ac))
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
            q, r, s = frames_from_two_splines(velocities, second_spline_points)
        end

        if(config.color==Color.RAINBOW)
            color_range = range(HSV(0,1,1), stop=HSV(360,1,1), length=size(spline_points, 2))
            rgb_colors = map(hsv-> convert(RGB, hsv), color_range)
            rainbow_colors = map(rgb->map(channel->Int(floor(channel*255)), (rgb.r, rgb.g, rgb.b)), rgb_colors)
        end

        # filter
        if(config.filter==Filter.TANGENTS)
            fixed_indices = [1, size(spline_points, 2)]
            if(config.backbone_type==BackboneType.CARTOON || config.color==Color.SECONDARY_STRUCTURE || config.color==Color.RESIDUE)
                prev_res_idx = sample_to_residue_indices[1]
                prev_ss = fragments(chain)[sample_to_residue_indices[1]].properties[:SS]
                for (i, resIdx) in enumerate(sample_to_residue_indices)
                    if(resIdx!=prev_res_idx)
                        residue = fragments(chain)[resIdx]
                        if(config.color==Color.RESIDUE)
                            push!(fixed_indices, i-1, i)
                        end
                        if((config.backbone_type==BackboneType.CARTOON || config.color==Color.SECONDARY_STRUCTURE)
                            && residue.properties[:SS]!=prev_ss)
                            push!(fixed_indices, i-1, i)
                            prev_ss = residue.properties[:SS]
                        end
                        prev_res_idx = resIdx
                    end
                end
            end
            sort!(fixed_indices)
            unique!(fixed_indices)

            remaining_indices = filter_points_threshold(spline_points, q, r, s, fixed_indices)
            spline_points = spline_points[:, remaining_indices]
            sample_to_residue_indices = sample_to_residue_indices[remaining_indices]
            q = q[:, remaining_indices]
            r = r[:, remaining_indices]
            s = s[:, remaining_indices]
            if(config.color==Color.RAINBOW)
                rainbow_colors = rainbow_colors[remaining_indices]
            end
        end


        log_info(types, "Type of spline points: ", typeof(spline_points))

        circles::Vector{PlainNonStdMesh{U}} = []

        framesA::Vector{ColoredMesh} = []
        framesB::Vector{ColoredMesh} = []
        push!(framesB, local_frame_mesh(spline_points[:, 1], q[:, 1], r[:, 1], s[:, 1]))

        # iterate and create vertices
        for current_index=axes(spline_points, 2) # start- and endcap bases are the first and last splinepoints #TODO offset wieder auf 2, wenn Endkappen wieder funktionieren

            # sanity check: frame should be orthogonal
            if(!approx_zero(dot(q[:, current_index], r[:, current_index])) || !approx_zero(dot(q[:, current_index], s[:, current_index])) || !approx_zero(dot(s[:, current_index], r[:, current_index])))
                log_info(frame_rotation, current_index, " wrong angles ", dot(q[:, current_index], r[:, current_index]), " ", dot(q[:, current_index], s[:, current_index]), " ", dot(s[:, current_index], r[:, current_index]), " # ", q[:, current_index], " ", r[:, current_index], " ", s[:, current_index])
            end

            # debug output
            frameA = local_frame_mesh(spline_points[:, current_index], q[:, current_index], r[:, current_index], s[:, current_index])
            push!(framesA, frameA)

            # generate cross-section
            if(config.backbone_type==BackboneType.BACKBONE)
                circle_points = create_circle_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, config.stick_radius)
            elseif(config.backbone_type==BackboneType.RIBBON)
                circle_points = create_ellipse_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, T(3)*config.stick_radius, config.stick_radius)
            elseif(config.backbone_type==BackboneType.CARTOON)
                residue = fragments(chain)[sample_to_residue_indices[current_index]]
                structure = residue.properties[:SS]
                if(structure==BiochemicalAlgorithms.SecondaryStructure.NONE)
                    circle_points = create_circle_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, config.stick_radius)
                elseif(structure==BiochemicalAlgorithms.SecondaryStructure.HELIX)
                    circle_points = create_ellipse_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, T(3)*config.stick_radius, T(1.5)*config.stick_radius)
                elseif(structure==BiochemicalAlgorithms.SecondaryStructure.SHEET) # TODO arrows
                    circle_points = create_rectangle_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], config.resolution, T(3)*config.stick_radius, T(0.5)*config.stick_radius)
                end
            end
            circle = PlainNonStdMesh(circle_points, Vector{Vector{Int}}(), Vector{NTuple{3, Int}}())
            #TODO extra geometry at SS boundaries?
            
            if(config.color==Color.UNIFORM)
                color!(circle, uniform_color)
            elseif(config.color==Color.CHAIN)
                color!(circle, chain_colors[chain_num])
            elseif(config.color==Color.RAINBOW)
                color!(circle, rainbow_colors[current_index])
            elseif(config.color==Color.SECONDARY_STRUCTURE)
                residue = fragments(chain)[sample_to_residue_indices[current_index]]
                color!(circle, structure_color_mapping[residue.properties[:SS]])
            elseif(config.color==Color.RESIDUE)
                residue = fragments(chain)[sample_to_residue_indices[current_index]]
                color!(circle, amino_acid_color_mapping[residue.name])
            end

            push!(circles, circle)
        end

        # m1 = circlesToSimpleMesh(circles)
        # export_mesh_to_ply("circles.ply", m1)

        spline_mesh = connect_circles_to_tube(circles)
        log_info(types, "Type of spline mesh: ", typeof(spline_mesh))

        # ----- debug export -----
        # cs = merge_multiple_meshes(circle_meshes)
        # export_mesh_representation_to_ply("circles.ply", Representation(cs))

        fs = reduce(merge, framesA)
        export_mesh_to_ply("framesA.ply", fs)

        fs = reduce(merge, framesB)
        export_mesh_to_ply("framesB.ply", fs)

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