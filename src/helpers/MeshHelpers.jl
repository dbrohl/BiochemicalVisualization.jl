export rotate_in_direction!

# Takes a 2D mesh and adds a z-coordinate of 0 to each point
function lift_into_3d(mesh::SimpleMesh)
    vertices_3d = map(v -> Point(v.coords..., 0), mesh.vertices)
    result = SimpleMesh(vertices_3d, mesh.topology)
    return result
end

function color!(mesh, new_color)
    mesh.colors = repeat([new_color], size(mesh.vertices, 2))
end

function translate!(mesh, vec::AbstractArray{T}) where T <: AbstractFloat 
    mesh.vertices .+= vec
end

# Assumes the object in the origin and rotates it, so that its new local z-axis points to direction. 
# Returns the new object. 
function rotate_in_direction!(mesh::Union{PlainMesh, PlainNonStdMesh}, direction::AbstractArray{T}) where T<:AbstractFloat
    if(direction[2]==0)
        rotationAxis = Vec(T(0),T(1),T(0))
    else
        rotationAxis = Vec(T(1), -direction[1]/direction[2], T(0)) * sign(direction[2]) # rotation axis is perpendicular to direction and lies in the xy-plane
    end
    rotationAngle = Meshes.∠(Vec{3,T}(direction...), Vec(T(0), T(0), T(1)))
    mesh.vertices = (mesh.vertices' * AngleAxis(rotationAngle, rotationAxis...))'
end

# Merges multiple mesh objects into one. The geometry (vertices/edges/faces) does not change. 
function merge_multiple_meshes(meshes::AbstractVector{PlainMesh{T}}) where {T}
    num_points = sum(map(m -> size(m.vertices, 2), meshes))
    num_connects = sum(map(m -> size(m.connections, 2), meshes))

    points = reduce(hcat, map(m->m.vertices, meshes))

    colors = Array{Tuple{Int64, Int64, Int64}}(undef, num_points)
    connects = Array{Integer, 2}(undef, 3, num_connects)

    point_offset = 0
    connect_offset = 0
    for m in meshes
        point_len = size(m.vertices, 2)
        connect_len = size(m.connections, 2)

        
        colors[point_offset+1 : point_offset+point_len] = m.colors
        

        connects[:, connect_offset+1 : connect_offset+connect_len] = m.connections
        connects[:, connect_offset+1 : connect_offset+connect_len] .+= point_offset

        point_offset += point_len
        connect_offset += connect_len
    end

    return PlainMesh(points, connects, colors)
end

# Adds faces between circles to create the surface of a tube. 
function connect_circles_to_tube(circles::AbstractVector{PlainNonStdMesh{T}}) where {T}

    # collect all points
    points = reduce(hcat, map(m->m.vertices, circles))
    colors = vcat(map(c -> c.colors, circles)...)

    resolution = size(circles[1].connections, 1)
    connections = Array{Integer, 2}(undef, 3, (length(circles)-1)*(2*resolution))
    
    offset = 0
    connection_i = 1
    prev_indices = nothing
    for c in circles

        # add connections between circles
        current_indices = (offset+1):(offset+nvertices(c))
        if(prev_indices !== nothing)
            @assert length(current_indices)==length(prev_indices)

            shift, flip = determine_offset(points[:, current_indices[1]], points[:, current_indices[2]], points[:, prev_indices])
            prev_indices = circshift(prev_indices, -shift)
            if(flip)

                log_info(circle_index_correction, "flip", shift, " ", flip, " ", prev_indices, " ", current_indices) # TODO 1c4k has a problem and flip is detected
                colors[current_indices] = repeat([(0, 150, 0)], length(current_indices))
                reverse!(prev_indices)
            end

            m1 = reduce(hcat, collect(t) for t in zip(current_indices, prev_indices, circshift(prev_indices, 1)))
            connections[:, connection_i:connection_i+resolution-1] = m1
            connection_i += resolution

            m2 = reduce(hcat, collect(t) for t in zip(current_indices, circshift(prev_indices, 1), circshift(current_indices, 1)))
            connections[:, connection_i:connection_i+resolution-1] = m2
            connection_i += resolution
        end

        prev_indices = current_indices
        offset += nvertices(c)
    end

    log_info(types, "Type of points in connect_tube: ", typeof(points))

    
    return PlainMesh(points, connections, colors)
end

function determine_offset(p1, p2, circle_points)
    # TODO

    distances1 = map(cp -> norm(p1 .- cp), eachcol(circle_points))
    nearest1 = argmin(distances1)

    distances2 = map(cp -> norm(p2 .- cp), eachcol(circle_points))
    distances2[nearest1] = max(distances2...)+1
    nearest2 = argmin(distances2)

    a = nearest2-nearest1
    while a<=0
        a += size(circle_points, 2)
    end

    
    if(a>size(circle_points, 2)/2)

        log_info(circle_index_correction, nearest1, " ", nearest2)
        return nearest1-1, true
    else 
        return nearest1-1, false
    end
end


# Generates a Hemisphere as a SimpleMesh. 
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
                # quad = connect((resolution*(i-2)+(j-1), resolution*(i-2)+j, resolution*(i-1)+j, resolution*(i-1)+j-1))
                # push!(connections, quad)
                a = connect((resolution*(i-2)+(j-1), resolution*(i-2)+j, resolution*(i-1)+j))
                b = connect((resolution*(i-2)+(j-1), resolution*(i-1)+j, resolution*(i-1)+j-1))
                push!(connections, a, b)

                if(j==resolution) # connect with beginning
                    # quad = connect((resolution*(i-2)+j, resolution*(i-2)+1, resolution*(i-1)+1,  resolution*(i-1)+j))
                    # push!(connections, quad)
                    a = connect((resolution*(i-2)+j, resolution*(i-2)+1, resolution*(i-1)+1))
                    b = connect((resolution*(i-2)+j, resolution*(i-1)+1,  resolution*(i-1)+j))
                    push!(connections, a, b)
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