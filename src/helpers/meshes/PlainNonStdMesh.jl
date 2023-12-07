"""
PlainNonStdMesh represents a 3D mesh that is described by its vertices, the connections between them and the vertexcolors. 
The faces are not necessarily triangles, but can include an arbitrary number of vertices. 
"""
mutable struct PlainNonStdMesh{T}
    vertices::Matrix{T} # 3 rows, n cols
    normals::Matrix{T}
    connections::Vector{Vector{Int}} # (supports m-gons instead of only triangles)
    colors::Vector{NTuple{3, Int}}
end


function PlainNonStdMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}, colors::AbstractVector{NTuple{3, Int}}) where {Dim, T, V, TP}
    if length(colors)!=length(mesh.vertices)
        throw(ErrorException("Different number of vertices and colors. "))
    end

    vertices = hcat([collect(v.coords.coords) for v in mesh.vertices]...)

    normals = similar(vertices)
    for i = 1:Meshes.nvertices(mesh)
        adjacent_face_normals = filter(f -> mesh.vertices[i] ∈ f, faces(mesh, 2))
        normals[:, i] = mean(adjacent_face_normals)
    end

    connects::Vector{Vector{Int}} = []

    for (i, f) in enumerate(Meshes.elements(Meshes.topology(mesh)))
        push!(connects, collect(f.indices))
    end

    PlainNonStdMesh{T}(vertices, normals, connects, colors)
end

import Base.==
function ==(a::PlainNonStdMesh{T}, b::PlainNonStdMesh{U}) where {T, U}
    return (T==U 
    && a.vertices==b.vertices 
    && a.normals==b.normals
    && a.connections==b.connections 
    && a.colors==b.colors)
end

function PlainNonStdMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}, color::NTuple{3, Int}) where {Dim, T, V, TP}
    PlainNonStdMesh(mesh, repeat([color], length(mesh.vertices)))
end

function PlainNonStdMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}) where {Dim, T, V, TP}
    PlainNonStdMesh(mesh, (0, 0, 255))
end

nvertices(m::PlainNonStdMesh) = size(m.vertices, 2)
