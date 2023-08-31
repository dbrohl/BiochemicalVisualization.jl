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
        # sphere_radius = 0.5
        # sphere_mesh = simplexify(Sphere{3, U}((0,0,0), sphere_radius))

        # spheres = map(a -> Translate(U.(a.r)...)((sphere_mesh)), c_alphas)
        # m1 = reduce(merge, spheres)
        # export_mesh_to_ply("c-alpha.ply", m1)




        # real backbone

        c_alpha_spline = CatmullRom(hcat(map(x->x.r, c_alphas)...))
        
        spline_points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, accelerations::AbstractMatrix{T} =  c_alpha_spline(vertices_per_unit)
        min_vel = minimum(norm, eachcol(velocities))
        max_vel = maximum(norm, eachcol(velocities))
        min_acc = minimum(norm, eachcol(accelerations))
        max_acc = maximum(norm, eachcol(accelerations))

        log_info(types, "Type of spline points: ", typeof(spline_points))
        #println("#splinepoints", size(spline_points, 2))

        # verts = [Point3(x...) for x in eachcol(spline_points)]
        # connects = [connect((a, b)) for (a, b) in zip(1:length(verts)-1, 2:length(verts))]
        # export_mesh_to_ply("spline.ply", SimpleMesh(verts, connects))

        circles::Vector{PlainNonStdMesh{U}} = []

        circle_meshes::Vector{PlainMesh{U}} = []
        frames::Vector{ColoredMesh} = []
        mesh_edges = Matrix{Int}(undef, 3, resolution)
        mesh_edges[1, :] = collect(1:resolution)
        mesh_edges[2, :] = circshift!(collect(1:resolution), -1)
        mesh_edges[3, :] = repeat([resolution+1], resolution)

        warncounter = 0
        count_bits::Vector{Bool} = []
        for i=2:size(spline_points, 2)-2 # start- and endcap bases are the first and last splinepoints
            tangent = spline_points[:, i+1] .- spline_points[:, i]
            tangent_length = norm(tangent)

            if(tangent_length < 10^-5) 
                # Duplicate points (in one location) generate a normal vector with length 0, which leads to a wrong orientation of the circle.
                # (It also would be duplicated and unnecessary)
                continue
            end

            # construct artificial frame
            if(abs(tangent[3]) < 10^-5)
                continue
            end
            normal = [1, 2, U(-tangent[1]-2*tangent[2])/tangent[3]]
            binormal = [U(0), 0, 0]
            binormal[1] = tangent[2]*normal[3] - tangent[3]*normal[2]
            binormal[2] = tangent[3]*normal[1] - tangent[1]*normal[3]
            binormal[3] = tangent[1]*normal[2] - tangent[2]*normal[1]

            tangent /= norm(tangent)
            normal /= norm(normal)
            binormal /= norm(binormal)

            circle_points = create_circle_in_local_frame(spline_points[:, i], normal, binormal, resolution, stick_radius)
            circle = PlainNonStdMesh(circle_points, Vector{Vector{Int}}(), Vector{NTuple{3, Int}}())
            circle_mesh = PlainMesh([circle_points spline_points[:, i]], mesh_edges, repeat([(255, 255, 255)], resolution+1))
            #frame = local_frame_mesh(spline_points[:, i], tangent, normal, binormal)
            relative_vel_size = (norm(velocities[:, i])-min_vel)/(max_vel-min_vel)
            relative_acc_size = (norm(accelerations[:, i])-min_acc)/(max_acc-min_acc)
            #println(relative_vel_size)
            frame = local_arrow_mesh(spline_points[:, i], velocities[:, i], (255-Int(round(relative_vel_size*255)), 255-Int(round(relative_vel_size*255)), 255))
            frame = merge(frame, local_arrow_mesh(spline_points[:, i], accelerations[:, i], (255, 255-Int(round(relative_acc_size*255)), 255-Int(round(relative_acc_size*255)))))
            distance_to_previous_center = min([norm(spline_points[:, i-1] - a) for a in eachcol(circle.vertices)]...)
            distance_to_current_center = min([norm(spline_points[:, i] - a) for a in eachcol(circle.vertices)]...)

            if(distance_to_previous_center< distance_to_current_center)
                warncounter = 20
                log_info(damaged_mesh, "Possible broken mesh, chain $chain_num spline_point $i")
            end

            if(i%100==0) 
                count_bits = []
                num::Int32 = i÷100
                while(num!=0)
                    push!(count_bits, convert(Bool, num & 0x00000001))
                    num = num >> 1
                end
                color!(circle, COUNT_COLOR_START)
            else


                if(warncounter>0)
                    color!(circle, BREAK_COLOR)
                    warncounter -= 1
                else
                    if(length(count_bits)!=0)
                        val = pop!(count_bits)
                        color!(circle, val ? COUNT_COLOR_1 : COUNT_COLOR_0)
                    else

                        color!(circle, in(i%100, [0,1]) ? COUNT_COLOR : chain_colors[chain_num])
                        # color!(circle, BACKGROUND_COLOR)
                    end
                end
            end
            #color!(circle, chain_colors[chain_num])

            push!(circles, circle)
            push!(circle_meshes, circle_mesh)
            push!(frames, frame)

        end

        # m1 = circlesToSimpleMesh(circles)
        # export_mesh_to_ply("circles.ply", m1)

        spline_mesh = connect_circles_to_tube(circles)

        log_info(types, "Type of spline mesh: ", typeof(spline_mesh))

        # ----- debug export -----
        cs = merge_multiple_meshes(circle_meshes)
        export_mesh_representation_to_ply("circles.ply", Representation(cs))

        fs = reduce(merge, frames)
        export_mesh_to_ply("frames.ply", fs)

        # ----- add hemispheres to both ends -----

    #     # position
    #     start_cap = deepcopy(cap_mesh)
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