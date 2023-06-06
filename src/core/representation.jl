export Representation

struct Representation{T <: Real}
    primitives::AbstractVector{GeometryBasics.GeometryPrimitive{3, T}}
    colors::AbstractVector{String}

    function Representation{T}(
            primitives=Vector{GeometryBasics.GeometryPrimitive{3, T}}(), 
            colors=Vector{String}()) where {T}
        new(primitives, colors)
    end
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