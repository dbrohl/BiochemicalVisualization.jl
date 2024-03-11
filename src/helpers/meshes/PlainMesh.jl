"""PlainMesh represents a 3D triangle mesh that is described by its vertices, the connections between them and the vertexcolors. """
mutable struct PlainMesh{T}
    vertices::Matrix{T} # 3 rows, n cols
    normals::Matrix{T}
    connections::Matrix{Int} # 3 rows, m cols, represents triangles
    colors::Vector{NTuple{3, Int}}
end



function PlainMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}, colors::AbstractVector{NTuple{3, Int}}) where {Dim, T, V, TP}
    if length(colors)!=length(mesh.vertices)
        throw(ErrorException("Different number of vertices and colors. "))
    end

    connects = Matrix{Int}(undef, 3, Meshes.nelements(Meshes.topology(mesh))) # converting the connections first will catch all non-triangles before they lead to problems with normals
    for (i, f) in enumerate(Meshes.elements(Meshes.topology(mesh)))
        if length(f.indices)!=3
            throw(ErrorException("mesh contained a face that contains $(length(f.indices)) instead of 3 vertices. "))
        end
        connects[:, i] = collect(f.indices)
    end

    vertices = hcat([collect(v.coords.coords) for v in mesh.vertices]...)
    normals = similar(vertices)
    for i = 1:Meshes.nvertices(mesh)
        adjacent_face_normals = map(f -> Meshes.normal(f), filter(f -> mesh.vertices[i] ∈ f, collect(Meshes.faces(mesh, 2))))
        normals[:, i] = mean(adjacent_face_normals)
    end

    PlainMesh{T}(vertices, normals, connects, colors)
end

function PlainMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}, color::NTuple{3, Int}) where {Dim, T, V, TP}
    PlainMesh(mesh, repeat([color], size(mesh.vertices, 1)))
end

function PlainMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}) where {Dim, T, V, TP}
    PlainMesh(mesh, (0, 0, 255))
end

import Base.==
function ==(a::PlainMesh{T}, b::PlainMesh{U}) where {T, U}
    return (T==U 
    && a.vertices==b.vertices 
    && a.normals==b.normals
    && a.connections==b.connections 
    && a.colors==b.colors)
end

nvertices(m::PlainMesh) = size(m.vertices, 2)
nconnections(m::PlainMesh) = size(m.connections, 2)