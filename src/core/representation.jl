export Representation

# A Representation can contain either a multiple collections of primitives and their colors or a colored mesh.
# The mesh has to contain triangle faces only.
struct Representation{T <: Real}
    primitives::Dict{String, AbstractVector{GeometryBasics.GeometryPrimitive{3,T}}}
    vertices::AbstractVector{T}
    normals::AbstractVector{T}
    connections::AbstractVector{Int}
    colors::Dict{String, AbstractVector{String}}

    function Representation{T}(;
        primitives=Dict{String, Vector{GeometryBasics.GeometryPrimitive{3,T}}}(),
        vertices=Vector{T}(),
        normals=Vector{T}(),
        connections=Vector{Int}(),
        colors=Dict{String, Vector{String}}()) where {T}
        
        new{T}(primitives, vertices, normals, connections, colors)
    end
end

function Representation(mesh::PlainMesh{T}) where T
    p = vec(mesh.vertices)
    n = vec(mesh.normals)
    c = vec(mesh.connections) .- 1

    colors = Vector{String}(undef, length(mesh.colors))
    for (i, c) in enumerate(mesh.colors)
        colors[i] = rgb_to_hex(c, prefix="#")
    end

    return Representation{T}(
        vertices=p, 
        normals = n,
        connections=c, 
        colors=Dict([("mesh", colors)]))

end

import Base.==
function ==(a::Representation{T}, b::Representation{U}) where {T, U}
    return (T==U 
    && a.primitives==b.primitives
    && a.vertices==b.vertices 
    && a.normals==b.normals
    && a.connections==b.connections 
    && a.colors==b.colors)
end



MsgPack.msgpack_type(::Type{GeometryBasics.Cylinder3{T}})      where {T} = MsgPack.StructType()
MsgPack.msgpack_type(::Type{GeometryBasics.Sphere{T}})         where {T} = MsgPack.StructType()
MsgPack.msgpack_type(::Type{Representation{T}}) where {T} = MsgPack.StructType()

# convenience constructors
GeometryBasics.Sphere(center::Vector3{T}, r::T) where {T<:Real} = GeometryBasics.Sphere{T}(GeometryBasics.Point3{T}(center...), r)
GeometryBasics.Cylinder(origin::Vector3{T}, extremity::Vector3{T}, radius::T) where {T<:Real} = 
    GeometryBasics.Cylinder(GeometryBasics.Point3{T}(origin...), GeometryBasics.Point3{T}(extremity...), radius)

center(s::GeometryBasics.Sphere)   = s.center
center(c::GeometryBasics.Cylinder) = c.origin + 0.5 * (c.extremity - c.origin)


