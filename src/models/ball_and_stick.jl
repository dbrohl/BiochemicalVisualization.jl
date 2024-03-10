export prepare_ball_and_stick_model

function prepare_ball_and_stick_model(
        chain::Chain{T}, 
        config::Union{Nothing, BallStickConfig{T}}=nothing; 
        fixed_color::Union{Nothing, NTuple{3, Int}}=nothing) where {T <: Real}

    if config===nothing
        config = BallStickConfig{T}(T(0.4), T(0.2), Color.ELEMENT)
    end

    sys = parent_system(chain)

    # geometry
    spheres = map(a -> GeometryBasics.Sphere(a.r, config.sphere_radius), atoms(chain))
    sticks = [(atom_by_idx(sys, b.a1), 
                atom_by_idx(sys, b.a2)) for b in bonds(chain)]
    midpoints = map(s -> (s[1].r + T(0.5)*(s[2].r - s[1].r)), sticks)
    cylinders = collect(Iterators.flatten(map(((s,m),) -> (
        GeometryBasics.Cylinder(s[1].r, m, config.stick_radius), 
        GeometryBasics.Cylinder(m, s[2].r, config.stick_radius)), zip(sticks, midpoints))))

    # colors
    if config.color==Color.UNIFORM || config.color==Color.CHAIN
        color_string = get_string_color(config.color, first(atoms(chain)), fixed_color)
        sphere_colors = repeat([color_string], length(spheres))
        cylinder_colors = repeat([color_string], length(cylinders))

    else

        # cache colors for all atoms 
        color_dict = Dict{Int, String}()
        for atom in eachatom(chain)
            color_dict[atom.idx] = get_string_color(config.color, atom, fixed_color)
        end

        sphere_colors = Vector{String}(undef, length(color_dict))
        for (i, atom_idx) in enumerate(atoms_df(chain).idx)
            sphere_colors[i] = color_dict[atom_idx]
        end

        cylinder_colors = collect(
                    Iterators.flatten(
                        map(
                            s -> (color_dict[s[1].idx], color_dict[s[2].idx]), 
                            sticks)
                        )
                    )
    end

    result = Representation{T}(
        primitives=Dict([("spheres", spheres), ("cylinders", cylinders)]), 
        colors=Dict([("spheres", sphere_colors), ("cylinders", cylinder_colors)]))
    
    
    return result
end

function prepare_ball_and_stick_model(
    ac::System{T}, config::Union{Nothing, BallStickConfig{T}}=nothing; 
    fixed_color::Union{Nothing, NTuple{3, Int}}=nothing) where {T<:Real}
    return handle_multichain_model(ac, config, fixed_color, prepare_ball_and_stick_model)
end