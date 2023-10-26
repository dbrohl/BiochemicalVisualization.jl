@enumx BackboneType begin
    BACKBONE
    RIBBON
    CARTOON
end
@enumx Color begin
    UNIFORM
    CHAIN
    RAINBOW
    SECONDARY_STRUCTURE
    RESIDUE
    ELEMENT
end # TODO charge? hydrophobicity?
@enumx Spline begin
    CATMULL_ROM
    CUBIC_B
end # TODO more options
@enumx ControlPoints begin
    C_ALPHA
    MID_POINTS
end
@enumx Frame begin
    RMF
    SECOND_SPLINE #implies Carson&Bugg control points
end # TODO hybrid?
@enumx Filter begin
    NONE
    ANGLE
end

struct BackboneConfig
    stick_radius
    resolution

    backbone_type
    color
    spline
    control_point_strategy
    frame
    filter
end