function prepare_ball_and_stick_model(
        ac::AbstractAtomContainer{T}; 
        sphere_radius=T(0.4), 
        stick_radius=T(0.2), resolution=30) where {T<:Real}

    start_time = now()

    spheres = map(a -> GeometryBasics.Sphere(a.r, sphere_radius), atoms(ac))
    sphere_colors = [element_color_web(e) for e in atoms_df(ac).element]

    sticks = [(atom_by_idx(ac, b.a1), 
            atom_by_idx(ac, b.a2)) for b in bonds(ac)]

    midpoints = map(s -> (s[1].r + T(0.5)*(s[2].r - s[1].r)), sticks)

    cylinders = collect(Iterators.flatten(map(((s,m),) -> (
        GeometryBasics.Cylinder(s[1].r, m, stick_radius), 
        GeometryBasics.Cylinder(m, s[2].r, stick_radius)), zip(sticks, midpoints))))
    cylinder_colors = collect(Iterators.flatten(
        map(s -> (element_color_web(s[1].element), element_color_web(s[2].element)), sticks)))

    result = Representation{T}(
        primitives=Dict([("spheres", spheres), ("cylinders", cylinders)]), 
        colors=Dict([("spheres", sphere_colors), ("cylinders", cylinder_colors)]))
    
    log_info(time_info, "Generated ball&stick representation in $((now()-start_time).value/1000) seconds. ")
    return result
end
