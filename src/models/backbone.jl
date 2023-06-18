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
    spline_points = map(v->Point(v), c_alpha_spline())

    circle_prototype = discretize(Sphere(Point(U(0),U(0)), stick_radius), RegularDiscretization(resolution))
    circle_mesh = lift_into_3d(circle_prototype)

    #TODO add first and last point as Balls

    circles = []
    for i=1:length(spline_points)-1
        normal = spline_points[i+1]-spline_points[i]
        if(sqrt(sum(normal.coords .^2)) ≈ 0) 
            # Duplicate points (in one location) generate a normal vector with length 0, which leads to a wrong orientation of the circle.
            # (It also would be duplicated and unnecessary)
            continue
        end
        circle = rotate_in_direction(circle_mesh, normal)
        circle = ColoredMesh(Translate(spline_points[i].coords...)(circle), (100, 0, 50))
        push!(circles, circle)
    end
    spline_mesh = merge_circles(circles)
    #spline_mesh = reduce(merge, circles)
    # TODO add end caps
    # TODO cylinder mesh gets distorded when too much rotation between circles
    # TODO spline connection

    return 0, ColoredMesh(spline_mesh, (255, 255, 255))


end
