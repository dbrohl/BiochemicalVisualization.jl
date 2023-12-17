"""
PlainNonStdMesh represents a 3D mesh that is described by its vertices, the connections between them and the vertexcolors. 
The faces are not necessarily triangles: they can include an arbitrary number of vertices. 
"""
mutable struct PlainNonStdMesh{T}
    vertices::Matrix{T} # 3 rows, n cols
    normals::Matrix{T}
    connections::Vector{Vector{Int}} # (supports m-gons instead of only triangles)
    colors::Vector{NTuple{3, Int}}
end

import Base.==
function ==(a::PlainNonStdMesh{T}, b::PlainNonStdMesh{U}) where {T, U}
    return (T==U 
    && a.vertices==b.vertices 
    && a.normals==b.normals
    && a.connections==b.connections 
    && a.colors==b.colors)
end

nvertices(m::PlainNonStdMesh) = size(m.vertices, 2)
