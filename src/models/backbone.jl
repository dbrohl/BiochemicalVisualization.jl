export rotation_test, prepare_backbone_model


FLIP_COLOR::NTuple{3, Int} = (200, 200, 200)
BREAK_COLOR::NTuple{3, Int} = (255, 0, 0)
COUNT_COLOR_START::NTuple{3, Int} = (255, 0, 255)
COUNT_COLOR_0::NTuple{3, Int} = (0, 0, 255)
COUNT_COLOR_1::NTuple{3, Int} = (0, 255, 255)
BACKGROUND_COLOR::NTuple{3, Int} = (60, 60, 60)

function rotation_test()
    stick_radius = 0.5
    resolution=30
    U = Float64

    circle_prototype = discretize(Sphere(Point(U(0),U(0)), stick_radius), RegularDiscretization(resolution))
    circle_mesh = PlainNonStdMesh(lift_into_3d(circle_prototype))

    all_points = [circle_mesh.vertices [0;0;0]]
    all_connects = [(1:resolution)'; circshift(1:resolution, 1)'; repeat([resolution+1], resolution)']
    all_colors = repeat([(100, 100, 100)], resolution+1)
    circle_std_mesh = PlainMesh(all_points, all_connects, all_colors)

    all_perm_meshes::Vector{PlainMesh{U}} = []
    permutations = [[1, 2, 3],
    [1, 3, 2],
    [2, 1, 3],
    [2, 3, 1],
    [3, 1, 2],
    [3, 2, 1]]

    for (height, perm) in enumerate(permutations)

        res = 16
        r = collect(range(0, 2*π, length=res+1))[1:end-1]
        directions = Matrix{U}(undef, 3, 0)
        colors = []

        for (i, z_angle) in enumerate(r)
            layer_radius = cos(z_angle)
            z = sin(z_angle)

            for (j, xy_angle) in enumerate(r)
                x = cos(xy_angle) * layer_radius
                y = sin(xy_angle) * layer_radius
                directions = [directions [x; y; z][perm]]
                push!(colors, convert.(Int, floor.(((x, y, z) .+ 1).*128)))
            end
        end
        size = 5

        finished_circles = [circle_std_mesh]

        for (i, (dir, color)) in enumerate(zip(eachcol(directions), colors))
            circle = deepcopy(circle_std_mesh)
            rotate_in_direction!(circle, [dir[1]; dir[2]; dir[3]])
            translate!(circle, [dir[1]; dir[2]; dir[3]] .* size)

            color!(circle, color)
            push!(finished_circles, circle)
        end
        temp = merge_multiple_meshes(finished_circles)
        translate!(temp, [0.0; 0; 12*(height-1)])
        push!(all_perm_meshes, temp)

    end

    export_mesh_representation_to_ply("rotation_test.ply", Representation(merge_multiple_meshes(all_perm_meshes)))
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
    ac::AbstractAtomContainer{T}; 
    stick_radius=T(0.2), resolution=30) where {T<:Real}

    start_time = now()
    U = Float64
    if(T <: AbstractFloat)
        U = T
    end

    log_info(types, "Types: ", T, " ", U)

    vertices_per_unit = resolution / 2*π*stick_radius
    # circle_prototype = discretize(Sphere(Point(U(0),U(0)), stick_radius), RegularDiscretization(resolution))
    # circle_mesh = PlainNonStdMesh(lift_into_3d(circle_prototype))
    cap_mesh = PlainMesh(Hemisphere(stick_radius, resolution, U))

    #log_info(types, "Types of proto meshes: ", typeof(circle_mesh), " ", typeof(cap_mesh))
    chain_meshes::Vector{PlainMesh{U}} = []
    chain_colors = map(c->map(channel->Int(channel*255), (c.r, c.g, c.b)), collect(distinguishable_colors(nchains(ac)+1))[2:end])
    for (chain_num, chain) in enumerate(eachchain(ac))

        c_alphas = filter(x -> x.element==Elements.C && x.name=="CA", atoms(chain))
        @assert length(c_alphas)>=2 # TODO was sonst?

        # c alpha export
        # sphere_radius = 0.05
        # sphere_mesh = simplexify(Sphere{3, U}((0,0,0), sphere_radius))

        # spheres = map(a -> Translate(U.(a.r)...)((sphere_mesh)), c_alphas)
        # spheres = map(s -> ColoredMesh(s, (0, 0, 0)), spheres)
        # m1 = reduce(merge, spheres)
        # export_mesh_to_ply("c-alpha.ply", m1)




        # real backbone
        c_alpha_spline = CatmullRom(hcat(map(x->x.r, c_alphas)...))
        #c_alpha_spline = CubicB(chain)
        
        spline_points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, accelerations::AbstractMatrix{T} = c_alpha_spline(vertices_per_unit)
        for i = axes(velocities, 2)
            velocities[:, i] = velocities[:, i] ./ norm(velocities[:, i])
        end

        # delete points with too small tangent vectors to avoid numerical problems
        mask = .!approx_zero.(map(t -> norm(t), eachcol(velocities)))
        spline_points = spline_points[:, mask]
        velocities = velocities[:, mask]

        filtered_indices = filter_points_threshold(spline_points, velocities)
        spline_points = spline_points[:, filtered_indices]
        velocities = velocities[:, filtered_indices]

        # sphere_radius = 0.2
        # sphere_mesh = discretize(Sphere{3, Float64}((0,0,0), sphere_radius), RegularDiscretization(6))
        # debug_mesh = reduce(BiochemicalVisualization.merge, map(a -> ColoredMesh(Translate(Float64.(a)...)((sphere_mesh)), (0, 255, 0)), eachcol(spline_points)))
        # export_mesh_to_ply("filtered_points.ply", debug_mesh)

        q::AbstractMatrix{T}, r::AbstractMatrix{T}, s::AbstractMatrix{T} = rmf(spline_points, velocities)
        # TODO short chains with <4 CA atoms
        
        #filtered_indices = no_filter(spline_points)

        log_info(types, "Type of spline points: ", typeof(spline_points))

        circles::Vector{PlainNonStdMesh{U}} = []

        framesA::Vector{ColoredMesh} = []
        framesB::Vector{ColoredMesh} = []
        push!(framesB, local_frame_mesh(spline_points[:, 1], q[:, 1], r[:, 1], s[:, 1]))
        mesh_edges = Matrix{Int}(undef, 3, resolution)
        mesh_edges[1, :] = collect(1:resolution)
        mesh_edges[2, :] = circshift!(collect(1:resolution), -1)
        mesh_edges[3, :] = repeat([resolution+1], resolution)

        warncounter = 0
        for i=1:size(spline_points, 2) # start- and endcap bases are the first and last splinepoints #TODO offset wieder auf 2, wenn Endkappen wieder funktionieren
            current_index = i

            # generate circle for backbone
            circle_points = create_circle_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], resolution, stick_radius)
            #circle_points = create_ellipse_in_local_frame(spline_points[:, current_index], r[:, current_index], s[:, current_index], resolution, stick_radius)
            circle = PlainNonStdMesh(circle_points, Vector{Vector{Int}}(), Vector{NTuple{3, Int}}())

            # debug output
            frameA = local_frame_mesh(spline_points[:, current_index], q[:, current_index], r[:, current_index], s[:, current_index])

            # sanity check: RMF should be orthogonal
            if(!approx_zero(dot(q[:, current_index], r[:, current_index])) || !approx_zero(dot(q[:, current_index], s[:, current_index])) || !approx_zero(dot(s[:, current_index], r[:, current_index])))
                log_info(misc, current_index, " wrong angles ", dot(q[:, current_index], r[:, current_index]), " ", dot(q[:, current_index], s[:, current_index]), " ", dot(s[:, current_index], r[:, current_index]), " # ", q[:, current_index], " ", r[:, current_index], " ", s[:, current_index])
            end


            # Rudimentary problem detection
            if(i>1)
                distance_to_previous_center = min([norm(spline_points[:, i-1] - a) for a in eachcol(circle.vertices)]...)
                distance_to_current_center = min([norm(spline_points[:, current_index] - a) for a in eachcol(circle.vertices)]...)

                if(distance_to_previous_center < distance_to_current_center)
                    warncounter = 20
                    log_info(damaged_mesh, "Possible broken mesh, chain $chain_num spline_point $current_index")
                end
            end

            if(warncounter>0)
                color!(circle, BREAK_COLOR)
                warncounter -= 1
            else
                color!(circle, chain_colors[chain_num])
            end
            #color!(circle, chain_colors[chain_num])

            push!(circles, circle)
            push!(framesA, frameA)
            #push!(framesB, frameB)

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