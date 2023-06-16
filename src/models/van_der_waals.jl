function prepare_van_der_waals_model(
    ac::AbstractAtomContainer{T}; resolution=30) where {T<:Real}
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

    sphere_mesh = discretize(Sphere{3, U}((0,0,0), U(1)), RegularDiscretization(resolution))
    spheres = map(atoms(ac)) do a
        scaling = max(a.radius, U(1.0))
        temp = Stretch(scaling)(sphere_mesh)
        temp = Translate(a.r.x, a.r.y, a.r.z)(temp)
        return temp        
    end

    sphere_colors = [element_color_rgb(e) for e in atoms_df(ac).element]
    sphere_meshes = map(zip(spheres, sphere_colors)) do (s,c)
        ColoredMesh(s, c)
    end
    meshesRepresentation = reduce(merge, sphere_meshes)

    return geometryBasicsRepresentation, meshesRepresentation
end