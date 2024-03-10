export Color, BallStickConfig, StickConfig, VdWConfig
@enumx Color begin
    UNIFORM
    CHAIN
    RAINBOW
    SECONDARY_STRUCTURE
    RESIDUE
    ELEMENT
end # TODO charge? hydrophobicity?

mutable struct BallStickConfig{T}
    sphere_radius::T
    stick_radius::T
    color::Color.T
end

mutable struct StickConfig{T}
    stick_radius::T
    color::Color.T
end

mutable struct VdWConfig
    color::Color.T
end