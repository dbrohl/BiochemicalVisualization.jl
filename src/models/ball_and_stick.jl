function prepare_ball_and_stick_model(
        ac::AbstractAtomContainer{T}; 
        sphere_radius=T(0.4), 
        stick_radius=T(0.2), resolution=30) where {T<:Real}

    generate_geometryBasics_representation(ac, sphere_radius=sphere_radius, stick_radius=stick_radius), generate_mesh(ac, sphere_radius=sphere_radius, stick_radius=stick_radius, resolution=resolution)
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
    stick_radius, resolution) where {T<:Real}

    #TODO Question: Why parameterized?

    # Meshes.jl uses only AbstractFloats. If T would be Int64, or T<:AbstractIrrational, this would lead to problems.
    U = Float64
    if(T <: AbstractFloat)
        U = T
    end

    # generate prototype meshes (conversion Primitive -> Mesh happens only once)
    sphere_mesh = discretize(Sphere{3, U}((0,0,0), sphere_radius), RegularDiscretization(resolution)) # default resolution in Meshes.jl is 50
    cylinder_mesh = discretize(CylinderSurface(stick_radius), RegularDiscretization(resolution, 2))

    # spheres for atoms
    spheres = map(atoms(ac)) do a
        Translate(a.r.x, a.r.y, a.r.z)(sphere_mesh)
    end

    sphere_colors = [element_color_rgb(e) for e in atoms_df(ac).element]
    final_sphere_meshes = map(zip(spheres, sphere_colors)) do (s,c)
        ColoredMesh(s, c)
    end


    # cylinders for bonds
    final_cylinder_meshes::Vector{ColoredMesh} = []
    for bond in bonds(ac)
        start_atom = atom_by_idx(ac, bond.a1)
        end_atom = atom_by_idx(ac, bond.a2)
        distance = sqrt(sum((end_atom.r-start_atom.r) .^ 2))
        midpoint = start_atom.r + 0.5*(end_atom.r-start_atom.r)
        
        if(distance>0)

            # from start_atom to midpoint
            first_cylinder = Stretch(U(1), U(1), U(0.5*distance))(cylinder_mesh)

            relative_end_pos = Vec(U.(midpoint - start_atom.r))
            if(relative_end_pos[2]==0)
                rotationAxis = Vec(U(0),U(1),U(0))
            else
                rotationAxis = Vec(1, -relative_end_pos[1]/relative_end_pos[2], 0) * sign(relative_end_pos[2]) # rotation axis is perpendicular to relative_end_pos and lies in the xy-plane
            end
            rotationAngle = Meshes.∠(relative_end_pos, Vec(U(0), U(0), U(1)))
            first_cylinder = Rotate(AngleAxis(rotationAngle, rotationAxis...))(first_cylinder)

            first_cylinder = Translate(start_atom.r.data)(first_cylinder)
            push!(final_cylinder_meshes, ColoredMesh(first_cylinder, element_color_rgb(start_atom.element)))




            # from midpoint to end_atom
            second_cylinder = Stretch(U(1), U(1), U(0.5*distance))(cylinder_mesh)

            relative_end_pos = Vec(U.(end_atom.r - midpoint))
            if(relative_end_pos[2]==0)
                rotationAxis = Vec(U(0),U(1),U(0))
            else
                rotationAxis = Vec(1, -relative_end_pos[1]/relative_end_pos[2], 0) * sign(relative_end_pos[2]) # rotation axis is perpendicular to relative_end_pos and lies in the xy-plane
            end
            rotationAngle = Meshes.∠(relative_end_pos, Vec(U(0), U(0), U(1)))
            second_cylinder = Rotate(AngleAxis(rotationAngle, rotationAxis...))(second_cylinder)

            second_cylinder = Translate(U.(midpoint)...)(second_cylinder)
            push!(final_cylinder_meshes, ColoredMesh(second_cylinder, element_color_rgb(end_atom.element)))
            
        end
    end

    # combine everything
    m1 = reduce(merge, final_sphere_meshes)
    m2 = reduce(merge, final_cylinder_meshes)

    merge(m1, m2)
    
end