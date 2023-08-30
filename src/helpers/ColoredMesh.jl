export ColoredMesh, merge
struct ColoredMesh{Dim,T,V<:AbstractVector{Point{Dim,T}},TP<:Topology} <: Mesh{Dim,T}
    vertices::V
    topology::TP
    colors::Vector{NTuple{3, Int}}
end

function ColoredMesh(mesh::SimpleMesh{Dim, T, V, TP}, colors::AbstractVector{NTuple{3, Int}}) where {Dim, T, V, TP}
    if(length(colors)!=length(mesh.vertices))
        return error("Different numbers of vertices and colors")
    end
    ColoredMesh{Dim, T, V, TP}(mesh.vertices, mesh.topology, colors)
end

function ColoredMesh(mesh::SimpleMesh{Dim, T, V, TP}, color::NTuple{3, Int}) where {Dim, T, V, TP}
    ColoredMesh{Dim, T, V, TP}(mesh.vertices, mesh.topology, repeat([color], length(mesh.vertices)))
end

function merge(m1::ColoredMesh, m2::ColoredMesh)
    merged_mesh = Base.merge(convert(SimpleMesh, m1), convert(SimpleMesh, m2))
    merged_colors = [m1.colors; m2.colors]
    ColoredMesh(merged_mesh, merged_colors)
end

function vertices(m::ColoredMesh)
    m.vertices
end

nvertices(m::ColoredMesh) = length(m.vertices)
Base.convert(::Type{<:SimpleMesh}, m::ColoredMesh) = SimpleMesh(m.vertices, m.topology)