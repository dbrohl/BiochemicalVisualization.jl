function prepare_backbone_model(
    ac::AbstractAtomContainer{T}; 
    stick_radius=T(0.2), resolution=30) where {T<:Real}

    U = Float64
    if(T <: AbstractFloat)
        U = T
    end

    c_alphas = filter(x -> x.element==Elements.C && x.name=="CA", atoms(ac))

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

    vertices_per_unit = resolution / 2*π*stick_radius
    spline_points = map(v->Point(v), c_alpha_spline(vertices_per_unit))

    circle_prototype = discretize(Sphere(Point(U(0),U(0)), stick_radius), RegularDiscretization(resolution))
    circle_mesh = lift_into_3d(circle_prototype)

    circles::Vector{ColoredMesh} = []
    for i=2:length(spline_points)-2
        normal = spline_points[i+1]-spline_points[i]
        normal_length = norm(normal)

        if(normal_length < 10^-5) 
            # Duplicate points (in one location) generate a normal vector with length 0, which leads to a wrong orientation of the circle.
            # (It also would be duplicated and unnecessary)
            continue
        end
        circle = rotate_in_direction(circle_mesh, normal)
        circle = ColoredMesh(Translate(spline_points[i].coords...)(circle), in(i%100, [0,1]) ? (255, 0,0) : (255, 255, 255))
        push!(circles, circle)
    end
    spline_mesh = merge_circles(circles)

    # add hemispheres to both ends

    cap_mesh = Hemisphere(stick_radius, resolution, U)
    start_cap = Translate(spline_points[1].coords...)(rotate_in_direction(cap_mesh, -(spline_points[2]-spline_points[1])))
    end_cap = Translate(spline_points[end-1].coords...)(rotate_in_direction(cap_mesh, -(spline_points[end-1]-spline_points[end])))

    start_cap = ColoredMesh(start_cap, (255, 255, 255))
    end_cap = ColoredMesh(end_cap, (255, 255, 255))

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

        # connect start_cap to tube
        if(i==1)
            a = connect((start_cap_base_indices[i], first_circle_indices[i], first_circle_indices[resolution], start_cap_base_indices[resolution]))
        else
            a = connect((start_cap_base_indices[i], first_circle_indices[i], first_circle_indices[i-1], start_cap_base_indices[i-1]))
        end
        push!(connections, a)

        # # connect end_cap to tube
        if(i==1)
            a = connect((end_cap_base_indices[i], last_circle_indices[i], last_circle_indices[resolution], end_cap_base_indices[resolution]))
        else
            a = connect((end_cap_base_indices[i], last_circle_indices[i], last_circle_indices[i - 1], end_cap_base_indices[i-1]))
        end
        push!(connections, a)
    end

    result = ColoredMesh(SimpleMesh(points, connections), colors)

    debug()
    
    return 0, result


end

function circle_example()
    # TODO cylinder mesh gets distorded when too much rotation between circles
    colors_r::Vector{Tuple{Int64, Int64, Int64}} = []
    colors_b ::Vector{Tuple{Int64, Int64, Int64}}= [] 
    for i=1:resolution
        v = 50 + (i-1)*155/resolution
        v = Int64(floor(v))
        push!(colors_r, (255, v, v))
        push!(colors_b, (v, v, 255))
    end

    circle_prototype = discretize(Sphere(Point(U(0),U(0)), 1), RegularDiscretization(resolution))
    circle_mesh = lift_into_3d(circle_prototype)


    c1 = ColoredMesh(deepcopy(circle_mesh), colors_r)


    c2 = Rotate(AngleAxis(Float32(-3*2*π/resolution), Float32(0.0), Float32(0.0), Float32(1.0)))(circle_mesh)
    c2 = Rotate(AngleAxis(Float32(π), Float32(1.0), Float32(0.0), Float32(0.0)))(c2)

    c2 = Translate(Float32(0), Float32(0), Float32(0.5))(c2)

    c2 = ColoredMesh(c2, colors_b)


    res = merge_circles([c1,c2])
end

function debug()
    c1 = [Point(-1.700227655690024, 0.7680342240904516, -1.9419433894881324), Point(-1.7555466050300998, 0.7403057401179092, -1.9975069500495923), Point(-1.8297949480037516, 0.7348318517124957, -2.0345688474827686), Point(-1.9101344834044225, 0.7525590710770069, -2.0467207157854204), Point(-1.982673779358539, 0.7904221892261989, -2.0318614018790755), Point(-2.0348701311823594, 0.8418743276517343, -1.992560217208101), Point(-2.0576983023007016, 0.8980189398090378, -1.935612693401249), Point(-2.047211103175144, 0.9491481165132516, -1.870865573135264), Point(-2.0052218736108522, 0.9864211637655724, -1.8095142188451083), Point(-1.9389909210471055, 1.0033932251442612, -1.7621668415599232), Point(-1.8599701848226915, 0.9971296749309522, -1.7370102379068937), Point(-1.781823063521697, 0.9687135405856557, -1.7383942080858306), Point(-1.7180619103145849, 0.9230582352636065, -1.7660794471810337), Point(-1.679711602597334, 0.8680579703163885, -1.8152789304203565), Point(-1.6734032687617013, 0.8132227938210722, -1.8774856219647629)] 
    c2 = [Point(-1.7094668829253186, 0.7980469086967135, -2.0093683087069474), Point(-1.7749430869609621, 0.782513943604564, -2.0582350452697797), Point(-1.852287042201011, 0.7925725293047841, -2.0870985442172723), Point(-1.9281252791399022, 0.826483473852971, -2.090968016518312), Point(-1.989344685251481, 0.8783832562219129, -2.06917440991757), Point(-2.025359868085701, 0.9392979308844372, -2.025486038180973), Point(-2.029943463415356, 0.9986947902049943, -1.9674570187567038), Point(-2.002302933318977, 1.046303583594511, -1.9051210939085617), Point(-1.9472175827853662, 1.0738923247926144, -1.8492567077595492), Point(-1.8742121599752048, 1.07669065963437, -1.809523325784792), Point(-1.7959099724831293, 1.054214729773159, -1.7927912149523701), Point(-1.7258501743057852, 1.0103508293718722, -1.8019535046855244), Point(-1.676146745292948, 0.9526834317676252, -1.8354259474411474), Point(-1.6553938564201667, 0.8911837472571487, -1.887420858786723), Point(-1.667179881831328, 0.8364856349793895, -1.9489478467530807)]

    con1 = connect((collect(1:length(c1))..., ))
    con2 = connect((collect(1:length(c2))..., ))

    m1 = SimpleMesh(c1, [con1])
    m2 = SimpleMesh(c2, [con2])

    colors_r::Vector{Tuple{Int64, Int64, Int64}} = []
    colors_b ::Vector{Tuple{Int64, Int64, Int64}}= [] 
    for i=1:length(c1)
        v = 50 + (i-1)*155/length(c1)
        v = Int64(floor(v))
        push!(colors_r, (255, v, v))
        push!(colors_b, (v, v, 255))
    end

    write_mesh_as_ply("flip_debug.ply", merge(ColoredMesh(m1, colors_r), ColoredMesh(m2, colors_b)))
end