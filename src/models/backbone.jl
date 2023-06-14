function prepare_backbone_model(
    ac::AbstractAtomContainer{T}; 
    stick_radius=T(0.2)) where {T<:Real}

    sphere_radius = 0.5
    U = Float64
    if(T <: AbstractFloat)
        U = T
    end

    # sphere_mesh = simplexify(Sphere{3, U}((0,0,0), sphere_radius))

    c_alphas = filter(x -> x.element==Elements.C && x.name=="CA", atoms(ac))

    # spheres = map(a -> deepcopy(sphere_mesh), c_alphas)
    # map(zip(spheres, c_alphas)) do (s, a)
    #     shift!(s, Vec{3, U}(a.r.x, a.r.y, a.r.z))
    # end

    # sphere_colors = [element_color_rgb(a.element) for a in c_alphas]
    # sphere_meshes = map(zip(spheres, sphere_colors)) do (s,c)
    #     ColoredMesh(s, c)
    # end

    c_alpha_spline = CatmullRom(map(x->x.r, c_alphas))

    spline_points = map(v->Point(v), c_alpha_spline())

    seg = Meshes.Chain(spline_points...)
    segment_mesh = ColoredMesh(simplexify(seg), (0,0,0))


    # m1 = reduce(merge, sphere_meshes)
    # m2 = reduce(merge, cylinder_meshes)
    # merge(m1, segment_mesh)
    return 0, segment_mesh


end
