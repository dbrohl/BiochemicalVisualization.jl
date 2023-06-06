function prepare_van_der_waals_model(
    ac::AbstractAtomContainer{T}) where {T<:Real}
    # BasicGeometry Version
    # todo: get vdw radii
    spheres = map(a -> GeometryBasics.Sphere(a.r, max(a.radius, T(1.0))), atoms(ac))
    sphere_colors = [element_color(e) for e in atoms_df(ac).element]

    # TODO Meshes Version

    Representation{T}(spheres, sphere_colors), 0
end