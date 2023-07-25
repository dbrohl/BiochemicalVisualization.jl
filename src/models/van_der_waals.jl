function prepare_van_der_waals_model(
    ac::AbstractAtomContainer{T}; resolution=30) where {T<:Real}

    # todo: get vdw radii
    start_time = now()
    spheres = map(a -> GeometryBasics.Sphere(a.r, max(a.radius, T(1.0))), atoms(ac))
    sphere_colors = [element_color_web(e) for e in atoms_df(ac).element]

    result = Representation{T}(primitives=Dict([("spheres", spheres)]), colors=Dict([("spheres", sphere_colors)]))
    log_info(time_info, "Generated vdw representation in $((now()-start_time).value/1000) seconds. ")
    
    return result 
end