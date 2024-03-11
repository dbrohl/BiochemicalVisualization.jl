export ColoredMesh, merge
struct ColoredMesh{Dim,T,V<:AbstractVector{Meshes.Point{Dim,T}},TP<:Meshes.Topology}
    vertices::V
    topology::TP
    colors::Vector{NTuple{3, Int}}
end

function ColoredMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}, colors::AbstractVector{NTuple{3, Int}}) where {Dim, T, V, TP}
    if(length(colors)!=length(mesh.vertices))
        return error("Different numbers of vertices and colors")
    end
    ColoredMesh{Dim, T, V, TP}(mesh.vertices, mesh.topology, colors)
end

function ColoredMesh(mesh::Meshes.SimpleMesh{Dim, T, V, TP}, color::NTuple{3, Int}) where {Dim, T, V, TP}
    ColoredMesh{Dim, T, V, TP}(mesh.vertices, mesh.topology, repeat([color], length(mesh.vertices)))
end

import Base.==
function ==(a::ColoredMesh{Dim,T,V,TP}, b::ColoredMesh{Dim2,T2,V2,TP2}) where {Dim, Dim2, T, T2, V, V2, TP, TP2}
    return (Dim==Dim2 && T==T2 && V==V2 && TP==TP2
    && a.vertices==b.vertices 
    && a.topology==b.topology 
    && a.colors==b.colors)
end

function merge(m1::ColoredMesh, m2::ColoredMesh)
    merged_mesh = Base.merge(convert(Meshes.SimpleMesh, m1), convert(Meshes.SimpleMesh, m2))
    merged_colors = [m1.colors; m2.colors]
    ColoredMesh(merged_mesh, merged_colors)
end

function vertices(m::ColoredMesh)
    m.vertices
end

nvertices(m::ColoredMesh) = length(m.vertices)
Base.convert(::Type{<:Meshes.SimpleMesh}, m::ColoredMesh) = Meshes.SimpleMesh(m.vertices, m.topology)