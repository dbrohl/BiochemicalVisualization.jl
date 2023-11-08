mutable struct PlainNonStdMesh{T}
    vertices::Matrix{T} # 3 rows, n cols
    connections::Vector{Vector{Int}} # (supports m-gons instead of only triangles)
    colors::Vector{NTuple{3, Int}}
end


function PlainNonStdMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}, colors::AbstractVector{NTuple{3, Int}}) where {Dim, T, V, TP}
    @assert length(colors)==length(mesh.vertices)

    vertices = hcat([collect(v.coords.coords) for v in mesh.vertices]...)

    connects::Vector{Vector{Int}} = []

    for (i, f) in enumerate(elements(topology(mesh)))
        push!(connects, collect(f.indices))
    end

    PlainNonStdMesh{T}(vertices, connects, colors)
end

function PlainNonStdMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}, color::NTuple{3, Int}) where {Dim, T, V, TP}
    PlainNonStdMesh(mesh, repeat([color], length(mesh.vertices)))
end

function PlainNonStdMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}) where {Dim, T, V, TP}
    PlainNonStdMesh(mesh, (0, 0, 255))
end

nvertices(m::PlainNonStdMesh) = size(m.vertices, 2)