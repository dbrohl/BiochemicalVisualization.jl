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
    LINEAR
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
mutable struct BackboneConfig{T}
    stick_radius::T
    resolution_along::T # the sampling frequency along the spline
    resolution_cross::Int # the number of vertices in one cross section
    backbone_type::BackboneType.T
    color::Color.T
    spline::Spline.T
    control_point_strategy::ControlPoints.T
    frame::Frame.T
    filter::Filter.T
end

import Base.==
function ==(a::BackboneConfig{T}, b::BackboneConfig{U}) where {T, U}
    return (T==U
    && a.stick_radius==b.stick_radius
    && a.resolution_along==b.resolution_along
    && a.resolution_cross==b.resolution_cross 
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
    resolution_along::Union{Real,  Nothing} # the sampling frequency along the spline
    resolution_cross::Union{Int,  Nothing} # the number of vertices in one cross section

    backbone_type::Union{BackboneType.T, Nothing}
    color::Union{Color.T, Nothing}
    spline::Union{Spline.T, Nothing}
    control_point_strategy::Union{ControlPoints.T, Nothing}
    frame::Union{Frame.T, Nothing}
    filter::Union{Filter.T, Nothing}
end

function PartialBackboneConfig(;    
    stick_radius = nothing, 
    resolution_along = nothing,
    resolution_cross = nothing, 
    backbone_type = nothing, 
    color = nothing, 
    spline = nothing, 
    control_point_strategy = nothing, 
    frame = nothing, 
    filter = nothing)
    return PartialBackboneConfig(stick_radius, resolution_along, resolution_cross, backbone_type, color, spline, control_point_strategy, frame, filter)
end

function ==(a::PartialBackboneConfig, b::PartialBackboneConfig)
    return (a.stick_radius==b.stick_radius
    && a.resolution_along==b.resolution_along 
    && a.resolution_cross==b.resolution_cross
    && a.backbone_type==b.backbone_type 
    && a.color==b.color
    && a.spline==b.spline
    && a.control_point_strategy==b.control_point_strategy
    && a.frame==b.frame
    && a.filter==b.filter)
end

function complete_config(partial::PartialBackboneConfig, template::BackboneConfig{T}) where {T}
    return BackboneConfig(
        partial.stick_radius===nothing ? template.stick_radius : T(partial.stick_radius),
        partial.resolution_along===nothing ? template.resolution_along : T(partial.resolution_along),
        partial.resolution_cross===nothing ? template.resolution_cross : partial.resolution_cross,

        partial.backbone_type===nothing ? template.backbone_type : partial.backbone_type,
        partial.color===nothing ? template.color : partial.color,
        partial.spline===nothing ? template.spline : partial.spline,
        partial.control_point_strategy===nothing ? template.control_point_strategy : partial.control_point_strategy,
        partial.frame===nothing ? template.frame : partial.frame,
        partial.filter===nothing ? template.filter : partial.filter,)
end