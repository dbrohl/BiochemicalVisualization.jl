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


# ----- misc ------
function lift_into_3d(mesh::SimpleMesh)
    vertices_3d = map(v -> Point(v.coords..., 0), mesh.vertices)
    result = SimpleMesh(vertices_3d, mesh.topology)
    return result
end

# Assumes the object in the origin and rotates it, so that its new local z-axis points to direction. 
# Returns the new object. 
function rotate_in_direction(mesh, direction::Vec{3, T}) where T<:AbstractFloat
    if(direction[2]==0)
        rotationAxis = Vec(T(0),T(1),T(0))
    else
        rotationAxis = Vec(1, -direction[1]/direction[2], 0) * sign(direction[2]) # rotation axis is perpendicular to direction and lies in the xy-plane
    end
    rotationAngle = Meshes.∠(direction, Vec(T(0), T(0), T(1)))
    result = Rotate(AngleAxis(rotationAngle, rotationAxis...))(mesh)
    return result
end

function merge_circles(circles)
    # collect all points
    points = vcat(map(c -> vertices(c), circles)...)

    connections::Vector{Connectivity} = []
    
    offset = 0
    prev_indices = nothing
    for c in circles
        # keep all existing circle-connections
        for primitive in elements(topology(c))
            a = connect(primitive.indices .+ offset)
            push!(connections, a)
        end

        # add connections between circles
        current_indices = (offset+1):(offset+nvertices(c))
        if(prev_indices !== nothing)
            @assert length(current_indices)==length(prev_indices)
            for i=1:length(current_indices)
                if(i==1)
                    a = connect((current_indices[i], prev_indices[i], prev_indices[end]))
                    b = connect((current_indices[i], current_indices[end], prev_indices[end]))
                else
                    a = connect((current_indices[i], prev_indices[i], prev_indices[i-1]))
                    b = connect((current_indices[i], current_indices[i-1], prev_indices[i-1]))
                end
                push!(connections, a, b)
            end
        end

        prev_indices = current_indices
        offset += nvertices(c)
    end

    return SimpleMesh(points, connections)
end