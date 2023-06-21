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

function merge_circles(circles::AbstractVector{T}) where T<: Union{SimpleMesh, ColoredMesh}
    # collect all points
    points::Vector{Point3} = vcat(map(c -> vertices(c), circles)...)
    colors = nothing
    if(T<:ColoredMesh)
        colors = vcat(map(c -> c.colors, circles)...)
    end

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

            shift, flip = determine_offset(points[current_indices[1]], points[current_indices[2]], points[prev_indices])
            prev_indices = circshift(prev_indices, -shift)
            if(flip)
                println("flip")
                reverse!(prev_indices)
            end
            for i in eachindex(current_indices)
                if(i==1)
                    a = connect((current_indices[i], prev_indices[i], prev_indices[end], current_indices[end]))
                else
                    a = connect((current_indices[i], prev_indices[i], prev_indices[i-1], current_indices[i-1]))
                end
                push!(connections, a)
            end
        end

        prev_indices = current_indices
        offset += nvertices(c)
    end

    if(colors===nothing)
        return SimpleMesh(points, connections)
    else
        return ColoredMesh(SimpleMesh(points, connections), colors)
    end
end

function determine_offset(p1, p2, circle_points)
    distances1 = map(cp -> norm(p1-cp), circle_points)
    nearest1 = argmin(distances1)

    distances2 = map(cp -> norm(p2-cp), circle_points)
    distances2[nearest1] = max(distances2...)+1
    nearest2 = argmin(distances2)

    a = nearest2-nearest1
    while a<=0
        a += length(circle_points)
    end

    
    if(a>length(circle_points)/2)

        println(nearest1, " ", nearest2)
        return nearest1-1, true
    else 
        return nearest1-1, false
    end
end



function Hemisphere(radius, resolution, T)

    vertices::Vector{Point{3,T}} = []
    connections::Vector{Connectivity} = []

    z_resolution = resolution ÷ 2
    xy_angles = collect(range(0, 2*π, resolution+1))[1: (end-1)]
    z_angles = collect(range(0, π/2, z_resolution+1))[1: (end-1)]

    
    for (i, z_angle) in enumerate(z_angles)
        layer_radius = cos(z_angle) * radius
        z = T(sin(z_angle) * radius)

        for (j, xy_angle) in enumerate(xy_angles)
            x = T(cos(xy_angle) * layer_radius)
            y = T(sin(xy_angle) * layer_radius)
            push!(vertices, Point(x, y, z))

            if(i>1 && j>1)
                quad = connect((resolution*(i-2)+(j-1), resolution*(i-2)+j, resolution*(i-1)+j, resolution*(i-1)+j-1))
                push!(connections, quad)

                if(j==resolution) # connect with beginning
                    quad = connect((resolution*(i-2)+j, resolution*(i-2)+1, resolution*(i-1)+1,  resolution*(i-1)+j))
                    push!(connections, quad)
                end

                if(i==z_resolution)
                    tri = connect((resolution*(i-1)+j-1,  resolution*(i-1)+j, resolution*z_resolution+1))
                    push!(connections, tri)

                    if(j==resolution)
                        tri = connect((resolution*(i-1)+j,  resolution*(i-1)+1, resolution*z_resolution+1))
                        push!(connections, tri)
                    end
                end
            end
        end
    end

    north = Point(0, 0, radius)
    push!(vertices, north)


    return SimpleMesh(vertices, connections)
end