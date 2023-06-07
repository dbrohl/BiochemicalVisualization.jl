function prepare_van_der_waals_model(
    ac::AbstractAtomContainer{T}) where {T<:Real}
    # BasicGeometry Version
    # todo: get vdw radii
    spheres = map(a -> GeometryBasics.Sphere(a.r, max(a.radius, T(1.0))), atoms(ac))
    sphere_colors = [element_color(e) for e in atoms_df(ac).element]

    geometryBasicsRepresentation = Representation{T}(spheres, sphere_colors)

    # Meshes.jl Version
    U = Float64
    if(T <: AbstractFloat)
        U = T
    end

    spheres = map(a -> Sphere{3, U}((a.r.x, a.r.y, a.r.z), max(a.radius, T(1.0))), atoms(ac))
    sphere_colors = [element_color_rgb(e) for e in atoms_df(ac).element]
    sphere_meshes = map(zip(spheres, sphere_colors)) do (s,c)
        ColoredMesh(simplexify(s), c)
    end
    meshesRepresentation = reduce(merge, sphere_meshes)

    return geometryBasicsRepresentation, meshesRepresentation
end