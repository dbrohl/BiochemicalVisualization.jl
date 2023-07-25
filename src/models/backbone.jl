export prepare_backbone_model
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
    circle_prototype = discretize(Sphere(Point(U(0),U(0)), stick_radius), RegularDiscretization(resolution))
    circle_mesh = PlainNonStdMesh(lift_into_3d(circle_prototype))
    cap_mesh = PlainMesh(Hemisphere(stick_radius, resolution, U))

    log_info(types, "Types of proto meshes: ", typeof(circle_mesh), " ", typeof(cap_mesh))
    chain_meshes::Vector{PlainMesh{U}} = []
    chain_colors = map(c->map(channel->Int(channel*255), (c.r, c.g, c.b)), collect(distinguishable_colors(nchains(ac)+1))[2:end])
    for (chain_num, chain) in enumerate(eachchain(ac))

        c_alphas = filter(x -> x.element==Elements.C && x.name=="CA", atoms(chain))
        @assert length(c_alphas)>=2 # TODO was sonst?

        # c alpha export
        # sphere_radius = 0.5
        # sphere_mesh = simplexify(Sphere{3, U}((0,0,0), sphere_radius))

        # spheres = map(a -> Translate(U.(a.r)...)((sphere_mesh)), c_alphas)

        # sphere_colors = [element_color_rgb(a.element) for a in c_alphas]
        # sphere_meshes = map(zip(spheres, sphere_colors)) do (s,c)
        #     ColoredMesh(s, c)
        # end
        # m1 = reduce(merge, sphere_meshes)
        # return 0, m1




        # real backbone

        c_alpha_spline = CatmullRom(hcat(map(x->x.r, c_alphas)...))

        
        spline_points::AbstractMatrix{T} =  c_alpha_spline(vertices_per_unit)
        log_info(types, "Type of spline points: ", typeof(spline_points))

        circles::Vector{PlainNonStdMesh{U}} = []
        for i=2:size(spline_points, 2)-2
            normal = spline_points[:, i+1] .- spline_points[:, i]
            normal_length = norm(normal)

            if(normal_length < 10^-5) 
                # Duplicate points (in one location) generate a normal vector with length 0, which leads to a wrong orientation of the circle.
                # (It also would be duplicated and unnecessary)
                continue
            end
            circle = deepcopy(circle_mesh)
            rotate_in_direction!(circle, normal)
            translate!(circle, spline_points[:, i])
            color!(circle, in(i%100, [0,1]) ? (255, 0,0) : chain_colors[chain_num])
            push!(circles, circle)
        end

        spline_mesh = connect_circles_to_tube(circles)

        log_info(types, "Type of spline mesh: ", typeof(spline_mesh))

        # ----- add hemispheres to both ends -----

        # position
        start_cap = deepcopy(cap_mesh)
        rotate_in_direction!(start_cap, -(spline_points[:, 2]-spline_points[:, 1]))
        translate!(start_cap, spline_points[:, 1])
        end_cap = deepcopy(cap_mesh)
        rotate_in_direction!(end_cap, -(spline_points[:, end-1]-spline_points[:, end]))
        translate!(end_cap, spline_points[:, end-1])

        # color!(start_cap, chain_colors[chain_num])
        # color!(end_cap, chain_colors[chain_num])
        color!(start_cap, (0, 255, 0))
        color!(end_cap, (0, 0, 255))

        log_info(types, "Typeof shifted cap: ", typeof(start_cap), " ", typeof(end_cap))

        # merge meshes

        points = [spline_mesh.vertices start_cap.vertices end_cap.vertices]
        colors = [spline_mesh.colors; start_cap.colors; end_cap.colors]

        connections::Matrix{Int} = Matrix(undef, 3, nconnections(spline_mesh) + nconnections(start_cap) + nconnections(end_cap) + 4 * resolution)
        
        position = 1
        offset = 0
        connections[:, position : position + nconnections(spline_mesh)-1] = spline_mesh.connections
        position += nconnections(spline_mesh)
        offset += nvertices(spline_mesh)

        connections[:, position : position + nconnections(start_cap)-1] = (start_cap.connections .+ offset)
        position += nconnections(start_cap)
        offset += nvertices(start_cap)

        connections[:, position : position + nconnections(end_cap)-1] = (end_cap.connections .+ offset)
        position += nconnections(end_cap)
        offset += nvertices(end_cap)

        # add new connections

        start_cap_base_indices = nvertices(spline_mesh)+1 : nvertices(spline_mesh)+resolution
        end_cap_base_indices = nvertices(spline_mesh)+nvertices(start_cap)+1 : nvertices(spline_mesh)+nvertices(start_cap)+resolution
        first_circle_indices = 1 : resolution
        last_circle_indices = nvertices(spline_mesh)-resolution+1 : nvertices(spline_mesh)

        shift, flip = determine_offset(points[:, start_cap_base_indices[1]], points[:, start_cap_base_indices[2]], points[:, first_circle_indices])
        first_circle_indices = circshift(first_circle_indices, -shift)
        if(flip)
            reverse!(first_circle_indices)
        end

        shift, flip = determine_offset(points[:, end_cap_base_indices[1]], points[:, end_cap_base_indices[2]], points[:, last_circle_indices])
        end_cap_base_indices = circshift(end_cap_base_indices, -shift)
        if(flip)
            reverse!(end_cap_base_indices)
        end


        connections[1, position:position+resolution-1] = start_cap_base_indices'
        connections[2, position:position+resolution-1] = first_circle_indices'
        connections[3, position:position+resolution-1] = circshift(first_circle_indices, 1)'
        position += resolution

        connections[1, position:position+resolution-1] = start_cap_base_indices'
        connections[2, position:position+resolution-1] = circshift(first_circle_indices, 1)'
        connections[3, position:position+resolution-1] = circshift(start_cap_base_indices, 1)'
        position += resolution

        connections[1, position:position+resolution-1] = end_cap_base_indices'
        connections[2, position:position+resolution-1] = last_circle_indices'
        connections[3, position:position+resolution-1] = circshift(last_circle_indices, 1)'
        position += resolution

        connections[1, position:position+resolution-1] = end_cap_base_indices'
        connections[2, position:position+resolution-1] = circshift(last_circle_indices, 1)'
        connections[3, position:position+resolution-1] = circshift(end_cap_base_indices, 1)'
        position += resolution
        
        push!(chain_meshes, PlainMesh(points, connections, colors))
        log_info(types, "Type of capped spline mesh: ", typeof(chain_meshes[end]))
    end
    temp = merge_multiple_meshes(chain_meshes)
    result = Representation(temp)
    log_info(types, "Type of result: ", typeof(result))

    log_info(time_info, "Generated backbone mesh in $((now()-start_time).value/1000) seconds. ($(length(result.vertices) ÷ 3) vertices)")

    return result

end