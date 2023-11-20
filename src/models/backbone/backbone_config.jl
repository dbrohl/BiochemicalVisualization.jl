export BackboneConfig, # TODO remove once the export of prepare_backbone_model is removed
PartialBackboneConfig, 
BackboneType, Color, Spline, ControlPoints, Frame, Filter

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

"Collection of parameters that determine how backbone-based vizualizations are created. "
mutable struct BackboneConfig
    stick_radius::Real
    resolution::Int

    backbone_type::BackboneType.T
    color::Color.T
    spline::Spline.T
    control_point_strategy::ControlPoints.T
    frame::Frame.T
    filter::Filter.T
end

import Base.==
function ==(a::BackboneConfig, b::BackboneConfig)
    return (a.stick_radius==b.stick_radius
    && a.resolution==b.resolution 
    && a.backbone_type==b.backbone_type 
    && a.color==b.color
    && a.spline==b.spline
    && a.control_point_strategy==b.control_point_strategy
    && a.frame==b.frame
    && a.filter==b.filter)
end

"Similar to BackboneConfig, but the datatypes include Nothing"
mutable struct PartialBackboneConfig
    stick_radius::Union{Real, Nothing}
    resolution::Union{Int, Nothing}

    backbone_type::Union{BackboneType.T, Nothing}
    color::Union{Color.T, Nothing}
    spline::Union{Spline.T, Nothing}
    control_point_strategy::Union{ControlPoints.T, Nothing}
    frame::Union{Frame.T, Nothing}
    filter::Union{Filter.T, Nothing}
end

function PartialBackboneConfig(;    
    stick_radius = nothing, 
    resolution = nothing, 
    backbone_type = nothing, 
    color = nothing, 
    spline = nothing, 
    control_point_strategy = nothing, 
    frame = nothing, 
    filter = nothing)
    return PartialBackboneConfig(stick_radius, resolution, backbone_type, color, spline, control_point_strategy, frame, filter)
end

function ==(a::PartialBackboneConfig, b::PartialBackboneConfig)
    return (a.stick_radius==b.stick_radius
    && a.resolution==b.resolution 
    && a.backbone_type==b.backbone_type 
    && a.color==b.color
    && a.spline==b.spline
    && a.control_point_strategy==b.control_point_strategy
    && a.frame==b.frame
    && a.filter==b.filter)
end

function complete_config(partial::PartialBackboneConfig, template::BackboneConfig)
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