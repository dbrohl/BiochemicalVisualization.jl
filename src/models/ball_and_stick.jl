function prepare_ball_and_stick_model(
        ac::AbstractAtomContainer{T}; 
        sphere_radius=T(0.4), 
        stick_radius=T(0.2)) where {T<:Real}

    generate_geometryBasics_representation(ac, sphere_radius=sphere_radius, stick_radius=stick_radius), generate_mesh(ac, sphere_radius=sphere_radius, stick_radius=stick_radius)
end

function generate_geometryBasics_representation(
    ac::AbstractAtomContainer{T}; 
    sphere_radius, 
    stick_radius) where {T<:Real}

    spheres = map(a -> GeometryBasics.Sphere(a.r, sphere_radius), atoms(ac))
    sphere_colors = [element_color(e) for e in atoms_df(ac).element]

    sticks = [(atom_by_idx(ac, b.a1), 
            atom_by_idx(ac, b.a2)) for b in bonds(ac)]

    midpoints = map(s -> (s[1].r + T(0.5)*(s[2].r - s[1].r)), sticks)

    cylinders = collect(Iterators.flatten(map(((s,m),) -> (
        GeometryBasics.Cylinder(s[1].r, m, stick_radius), 
        GeometryBasics.Cylinder(m, s[2].r, stick_radius)), zip(sticks, midpoints))))
    cylinder_colors = collect(Iterators.flatten(
        map(s -> (element_color(s[1].element), element_color(s[2].element)), sticks)))

    Representation{T}(vcat(spheres, cylinders), vcat(sphere_colors, cylinder_colors))

end

function generate_mesh(
    ac::AbstractAtomContainer{T}; 
    sphere_radius, 
    stick_radius) where {T<:Real}

    #TODO Question: Why parameterized?
    U = Float64
    if(T <: AbstractFloat)
        U = T
    end

    spheres = map(a -> Sphere{3, U}((a.r.x, a.r.y, a.r.z), sphere_radius), atoms(ac))
    sphere_colors = [element_color_rgb(e) for e in atoms_df(ac).element]
    sphere_meshes = map(zip(spheres, sphere_colors)) do (s,c)
        ColoredMesh(simplexify(s), c)
    end

    sticks = [(atom_by_idx(ac, b.a1), 
            atom_by_idx(ac, b.a2)) for b in bonds(ac)]

    midpoints = map(s -> (s[1].r + (0.5)*(s[2].r - s[1].r)), sticks)

    cylinders = collect(Iterators.flatten(map(((s,m),) -> (
        CylinderSurface(stick_radius, Segment(Point{3,U}(s[1].r), Point{3,U}(m))), 
        CylinderSurface(stick_radius, Segment(Point{3,U}(m), Point{3,U}(s[2].r)))), zip(sticks, midpoints))))

    cylinder_colors = collect(Iterators.flatten(
        map(s -> (element_color_rgb(s[1].element), element_color_rgb(s[2].element)), sticks)))
    
    cylinder_meshes = map(zip(cylinders, cylinder_colors)) do (cy,col)
        ColoredMesh(simplexify(cy), col)
    end

    m1 = reduce(merge, sphere_meshes)
    m2 = reduce(merge, cylinder_meshes)


    merge(m1, m2)

end