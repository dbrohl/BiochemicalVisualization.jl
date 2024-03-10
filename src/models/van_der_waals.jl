export prepare_van_der_waals_model

function prepare_van_der_waals_model(
    ac::System{T}, config::Union{Nothing, VdWConfig}=nothing; 
    fixed_color::Union{Nothing, NTuple{3, Int}}=nothing) where {T<:Real}
    return handle_multichain_model(ac, config, fixed_color, prepare_van_der_waals_model)
end

function prepare_van_der_waals_model(
    chain::Chain{T}, 
    config::Union{Nothing, VdWConfig}=nothing; 
    fixed_color::Union{Nothing, NTuple{3, Int}}=nothing) where {T <: Real}

    spheres = map(a -> GeometryBasics.Sphere(a.r, a.radius==0 ? T(1.0) : a.radius), atoms(chain))
    sphere_colors = [get_string_color(config.color, atom, fixed_color) for atom in eachatom(chain)]

    result = Representation{T}(primitives=Dict([("spheres", spheres)]), colors=Dict([("spheres", sphere_colors)]))
    
    return result 
end