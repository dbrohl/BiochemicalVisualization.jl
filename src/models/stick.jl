export prepare_stick_model

function prepare_stick_model(
    ac::Union{System{T}, Chain{T}}, 
    config::Union{Nothing, StickConfig{T}}=nothing; 
    fixed_color::Union{Nothing, NTuple{3, Int}}=nothing) where {T <: Real}

    return prepare_ball_and_stick_model(ac, BallStickConfig{T}(config.stick_radius, config.stick_radius, config.color), fixed_color=fixed_color)
end