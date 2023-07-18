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
    circle_mesh = lift_into_3d(circle_prototype)
    cap_mesh = Hemisphere(stick_radius, resolution, U)

    log_info(types, "Types of proto meshes: ", typeof(circle_mesh), " ", typeof(cap_mesh))
    chain_meshes::Vector{ColoredMesh{3, U}} = []
    chain_colors = map(c->map(channel->Int64(channel*255), (c.r, c.g, c.b)), collect(distinguishable_colors(nchains(ac)+1))[2:end])
    for (chain_num, chain) in enumerate(eachchain(ac))

        c_alphas = filter(x -> x.element==Elements.C && x.name=="CA", atoms(chain))
        @assert length(c_alphas)>=2

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

        c_alpha_spline = CatmullRom(map(x->Point(x.r), c_alphas))

        
        spline_points = map(v->Point(U.(v)), c_alpha_spline(vertices_per_unit))
        log_info(types, "Type of spline points: ", typeof(spline_points))

        circles::Vector{ColoredMesh{3, U}} = []
        for i=2:length(spline_points)-2
            normal = spline_points[i+1]-spline_points[i]
            normal_length = norm(normal)

            if(normal_length < 10^-5) 
                # Duplicate points (in one location) generate a normal vector with length 0, which leads to a wrong orientation of the circle.
                # (It also would be duplicated and unnecessary)
                continue
            end
            circle = rotate_in_direction(circle_mesh, normal)
            circle = ColoredMesh(Translate(spline_points[i].coords...)(circle), in(i%100, [0,1]) ? (255, 0,0) : chain_colors[chain_num])
            push!(circles, circle)
        end

        spline_mesh = connect_circles_to_tube(circles)

        log_info(types, "Type of spline mesh: ", typeof(spline_mesh))

        # add hemispheres to both ends

        
        start_cap = Translate(spline_points[1].coords...)(rotate_in_direction(cap_mesh, -(spline_points[2]-spline_points[1])))
        end_cap = Translate(spline_points[end-1].coords...)(rotate_in_direction(cap_mesh, -(spline_points[end-1]-spline_points[end])))

        start_cap = ColoredMesh(start_cap, chain_colors[chain_num])
        end_cap = ColoredMesh(end_cap, chain_colors[chain_num])

        log_info(types, "Typeof shifted cap: ", typeof(start_cap), " ", typeof(end_cap))

        start_cap_base_indices = nvertices(spline_mesh)+1 : nvertices(spline_mesh)+resolution
        end_cap_base_indices = nvertices(spline_mesh)+nvertices(start_cap)+1 : nvertices(spline_mesh)+nvertices(start_cap)+resolution
        first_circle_indices = 1 : resolution
        last_circle_indices = nvertices(spline_mesh)-resolution+1 : nvertices(spline_mesh)

        points = [vertices(spline_mesh); vertices(start_cap); vertices(end_cap)]
        colors = [spline_mesh.colors; start_cap.colors; end_cap.colors]

        connections::Vector{Connectivity} = []
        for primitive in elements(topology(spline_mesh))
            a = connect(primitive.indices)
            push!(connections, a)
        end
        for primitive in elements(topology(start_cap))
            a = connect(primitive.indices .+ nvertices(spline_mesh))
            push!(connections, a)
        end
        for primitive in elements(topology(end_cap))
            a = connect(primitive.indices .+ (nvertices(spline_mesh)+nvertices(start_cap)))
            push!(connections, a)
        end

        shift, flip = determine_offset(points[start_cap_base_indices[1]], points[start_cap_base_indices[2]], points[first_circle_indices])
        first_circle_indices = circshift(first_circle_indices, -shift)
        if(flip)
            reverse!(first_circle_indices)
        end

        shift, flip = determine_offset(points[end_cap_base_indices[1]], points[end_cap_base_indices[2]], points[last_circle_indices])
        end_cap_base_indices = circshift(end_cap_base_indices, -shift)
        if(flip)
            reverse!(end_cap_base_indices)
        end

        for i in 1:resolution # the circles and cap-bases must have the same number of vertices to connect them 1:1

            if(i==1)
                a = connect((start_cap_base_indices[i], first_circle_indices[i], first_circle_indices[resolution]))
                b = connect((start_cap_base_indices[i], first_circle_indices[resolution], start_cap_base_indices[resolution]))
            else
                a = connect((start_cap_base_indices[i], first_circle_indices[i], first_circle_indices[i-1]))
                b = connect((start_cap_base_indices[i], first_circle_indices[i-1], start_cap_base_indices[i-1]))
            end
            push!(connections, a, b)

            if(i==1)
                a = connect((end_cap_base_indices[i], last_circle_indices[i], last_circle_indices[resolution]))
                b = connect((end_cap_base_indices[i], last_circle_indices[resolution], end_cap_base_indices[resolution]))
            else
                a = connect((end_cap_base_indices[i], last_circle_indices[i], last_circle_indices[i - 1]))
                b = connect((end_cap_base_indices[i], last_circle_indices[i - 1], end_cap_base_indices[i-1]))
            end
            push!(connections, a, b)
        end
        
        push!(chain_meshes, ColoredMesh(SimpleMesh(points, connections), colors))
        log_info(types, "Type of capped spline mesh: ", typeof(chain_meshes[end]))
    end
    
    result = Representation(merge_multiple_meshes(chain_meshes))
    log_info(types, "Type of result: ", typeof(result))

    log_info(time_info, "Generated backbone mesh in $((now()-start_time).value/1000) seconds. ($(length(result.vertices) ÷ 3) vertices)")

    return result

end