export BackboneConfig, BackboneType, Color, Spline, ControlPoints, Frame, Filter

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

function BackboneConfig(;    
    stick_radius = nothing, 
    resolution = nothing, 
    backbone_type = nothing, 
    color = nothing, 
    spline = nothing, 
    control_point_strategy = nothing, 
    frame = nothing, 
    filter = nothing)
    return BackboneConfig(stick_radius, resolution, backbone_type, color, spline, control_point_strategy, frame, filter)
end

function complete_config(partial, template)
    return BackboneConfig(
        partial.stick_radius===nothing ? template.stick_radius : partial.stick_radius,
        partial.resolution===nothing ? template.resolution : partial.resolution,

        partial.backbone_type===nothing ? template.backbone_type : partial.backbone_type,
        partial.color===nothing ? template.color : partial.color,
        partial.spline===nothing ? template.spline : partial.spline,
        partial.control_point_strategy===nothing ? template.control_point_strategy : partial.control_point_strategy,
        partial.frame===nothing ? template.frame : partial.frame,
        partial.filter===nothing ? template.filter : partial.filter,)
end