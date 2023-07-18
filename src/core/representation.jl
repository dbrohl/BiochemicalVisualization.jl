export Representation

# A Representation can contain either a multiple collections of primitives and their colors or a colored mesh.
# The mesh has to contain triangle faces only.
struct Representation{T <: Real}
    primitives::Dict{String, AbstractVector{GeometryBasics.GeometryPrimitive{3,T}}}
    vertices::AbstractVector{T}
    connections::AbstractVector{Int}
    colors::Dict{String, AbstractVector{String}}

    function Representation{T}(;
        primitives=Dict{String, Vector{GeometryBasics.GeometryPrimitive{3,T}}}(),
        vertices=Vector{T}(),
        connections=Vector{Int}(),
        colors=Dict{String, Vector{String}}()) where {T}
        
        new{T}(primitives, vertices, connections, colors)
    end
end

function Representation(mesh::Union{SimpleMesh{Dim, T}, ColoredMesh{Dim, T}}) where {Dim, T}

    num_point_coords = length(mesh.vertices) * 3
    points = Array{T}(undef, num_point_coords)
    colors = Array{String}(undef, length(mesh.vertices))

    for i = 1:length(mesh.vertices)
        points[3*(i-1)+1] = mesh.vertices[i].coords[1]
        points[3*(i-1)+2] = mesh.vertices[i].coords[2]
        points[3*(i-1)+3] = mesh.vertices[i].coords[3]

        if(typeof(mesh)<:ColoredMesh)
            colors[i] = "#"*hex(RGB((mesh.colors[i] ./ 255)...))
        else
            colors[i] = "#0000ff"
        end
        
    end

    num_connects = nelements(topology(mesh)) * 3
    connections = Array{Int64}(undef, num_connects)  
    a = 0
    for f in elements(topology(mesh))
        @assert length(f.indices)==3
        connections[a+1:a+3] = [convert.(Int64, (f.indices .-1))...]
        a+=3
    end

    return Representation{T}(vertices=points, connections=connections, colors=Dict([("mesh", colors)]))
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


