mutable struct PlainNonStdMesh{T}
    vertices::AbstractMatrix{T} # 3 rows, n cols
    connections::Vector{Vector{Int}} # (supports m-gons instead of only triangles)
    colors::Vector{NTuple{3, Int}}
end


function PlainNonStdMesh(mesh::SimpleMesh{Dim, T, V, TP}, colors::AbstractVector{NTuple{3, Int}}) where {Dim, T, V, TP}
    @assert length(colors)==length(mesh.vertices)

    vertices = hcat([collect(v.coords.coords) for v in mesh.vertices]...)

    connects::Vector{Vector{Int}} = []

    for (i, f) in enumerate(elements(topology(mesh)))
        push!(connects, collect(f.indices))
    end

    PlainNonStdMesh{T}(vertices, connects, colors)
end

function PlainNonStdMesh(mesh::SimpleMesh{Dim, T, V, TP}, color::NTuple{3, Int}) where {Dim, T, V, TP}
    PlainNonStdMesh(mesh, repeat([color], length(mesh.vertices)))
end

function PlainNonStdMesh(mesh::SimpleMesh{Dim, T, V, TP}) where {Dim, T, V, TP}
    PlainNonStdMesh(mesh, (0, 0, 255))
end

function PlainNonStdMesh(mesh::ColoredMesh{Dim, T, V, TP}) where {Dim, T, V, TP}
    PlainNonStdMesh(convert(SimpleMesh, mesh), mesh.colors)
end


# function merge(m1::ColoredMesh, m2::ColoredMesh)
#     merged_mesh = Base.merge(convert(SimpleMesh, m1), convert(SimpleMesh, m2))
#     merged_colors = [m1.colors; m2.colors]
#     ColoredMesh(merged_mesh, merged_colors)
# end

# function vertices(m::PlainNonStdMesh)
#     m.vertices
# end

nvertices(m::PlainNonStdMesh) = size(m.vertices, 2)
# Base.convert(::Type{<:SimpleMesh}, m::ColoredMesh) = SimpleMesh(m.vertices, m.topology)