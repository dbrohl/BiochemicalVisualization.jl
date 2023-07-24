mutable struct PlainMesh{T}
    vertices::AbstractMatrix{T} # 3 rows, n cols
    connections::AbstractMatrix{Integer} # 3 rows, m cols, represents triangles
    colors::Vector{Tuple{Int64, Int64, Int64}}
end



function PlainMesh(mesh::SimpleMesh{Dim, T, V, TP}, colors::AbstractVector{NTuple{3, W}}) where {Dim, T, V, TP, W <: Integer}
    @assert length(colors)==length(mesh.vertices)

    vertices = reduce(hcat, [collect(v.coords.coords) for v in mesh.vertices])

    connects = Array{Int64, 2}(undef, 3, nelements(topology(mesh)))

    for (i, f) in enumerate(elements(topology(mesh)))
        @assert length(f.indices)==3
        connects[:, i] = collect(f.indices)
    end

    PlainMesh{T}(vertices, connects, colors)
end

function PlainMesh(mesh::SimpleMesh{Dim, T, V, TP}, color::NTuple{3, W}) where {Dim, T, V, TP, W <: Integer}
    PlainMesh(mesh, repeat([color], size(mesh.vertices, 1)))
end

function PlainMesh(mesh::SimpleMesh{Dim, T, V, TP}) where {Dim, T, V, TP, W <: Integer}
    PlainMesh(mesh, (0, 0, 255))
end

function PlainMesh(mesh::ColoredMesh{Dim, T, V, TP}) where {Dim, T, V, TP, W <: Integer}
    PlainMesh(convert(SimpleMesh, mesh), mesh.colors)
end



# function merge(m1::ColoredMesh, m2::ColoredMesh)
#     merged_mesh = Base.merge(convert(SimpleMesh, m1), convert(SimpleMesh, m2))
#     merged_colors = [m1.colors; m2.colors]
#     ColoredMesh(merged_mesh, merged_colors)
# end

# function vertices(m::ColoredMesh)
#     m.vertices
# end

nvertices(m::ColoredMesh) = size(m.vertices, 2)
# Base.convert(::Type{<:SimpleMesh}, m::ColoredMesh) = SimpleMesh(m.vertices, m.topology)