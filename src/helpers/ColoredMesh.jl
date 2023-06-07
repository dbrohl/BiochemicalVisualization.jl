struct ColoredMesh{Dim,T,V<:AbstractVector{Point{Dim,T}},TP<:Topology} <: Mesh{Dim,T} # TODO Why is it impossible to extend SimpleMesh?
    vertices::V
    topology::TP
    colors #::AbstractVector{NTuple{3, Integer}}
end

function ColoredMesh(mesh::SimpleMesh, colors::AbstractVector{NTuple{3, T}}) where T<:Integer
    if(length(colors)!=length(mesh.vertices))
        return error("Different numbers of vertices and colors")
    end
    ColoredMesh(mesh.vertices, mesh.topology, colors)
end

function ColoredMesh(mesh::SimpleMesh, color::NTuple{3, T}) where T<:Integer
    ColoredMesh(mesh.vertices, mesh.topology, repeat([color], length(mesh.vertices)))
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