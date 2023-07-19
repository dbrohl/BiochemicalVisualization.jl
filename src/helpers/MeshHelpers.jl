# Takes a 2D mesh and adds a z-coordinate of 0 to each point
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

# Merges multiple mesh objects into one. The geometry (vertices/edges/faces) does not change. 
function merge_multiple_meshes(meshes::AbstractVector{X}) where {Dim, T, X<:Union{SimpleMesh{Dim, T}, ColoredMesh{Dim, T}}}
    num_points = sum(map(m -> length(m.vertices), meshes))
    num_connects = sum(map(m -> nelements(topology(m)), meshes))

    points = Array{Point{3, T}}(undef, num_points)
    if(X<:ColoredMesh)
        colors = Array{Tuple{Int64, Int64, Int64}}(undef, num_points)
    end
    connects = Array{Connectivity}(undef, num_connects)

    point_offset = 0
    connect_offset = 0
    for m in meshes
        point_len = length(m.vertices)

        for i=1:point_len
            points[point_offset+i] = m.vertices[i]
        end

        if(X<:ColoredMesh)
            colors[point_offset+1 : point_offset+point_len] = m.colors
        end

        for primitive in elements(topology(m))
            connects[connect_offset+1] = connect(primitive.indices .+ point_offset)
            connect_offset += 1
        end

        point_offset += point_len
    end

    if(X<:ColoredMesh)
        return ColoredMesh(SimpleMesh(points, connects), colors)
    else
        return SimpleMesh(points, connects)
    end
end

# Adds faces between circles to create the surface of a tube. 
function connect_circles_to_tube(circles::AbstractVector{X}) where {Dim, T, X<: Union{ColoredMesh{Dim, T}, SimpleMesh{Dim, T}}}

    # collect all points
    points::Vector{Point{3, T}} = vcat(map(c -> vertices(c), circles)...)
    colors = nothing
    if(X<:ColoredMesh)
        colors = vcat(map(c -> c.colors, circles)...)
    end

    connections::Vector{Connectivity} = []
    
    offset = 0
    prev_indices = nothing
    for c in circles

        # add connections between circles
        current_indices = (offset+1):(offset+nvertices(c))
        if(prev_indices !== nothing)
            @assert length(current_indices)==length(prev_indices)

            shift, flip = determine_offset(points[current_indices[1]], points[current_indices[2]], points[prev_indices])
            prev_indices = circshift(prev_indices, -shift)
            if(flip)

                log_info(circle_index_correction, "flip", shift, " ", flip, " ", prev_indices, " ", current_indices) # TODO 1c4k has a problem and flip is detected
                # if(X<:ColoredMesh)
                #     colors[current_indices] = repeat([(0, 0, 255)], length(current_indices))
                # end
                reverse!(prev_indices)
            end
            for i in eachindex(current_indices)
                if(i==1)
                    a = connect((current_indices[i], prev_indices[i], prev_indices[end]))
                    b = connect((current_indices[i], prev_indices[end], current_indices[end]))
                else
                    a = connect((current_indices[i], prev_indices[i], prev_indices[i-1]))
                    b = connect((current_indices[i], prev_indices[i-1], current_indices[i-1]))
                end
                push!(connections, a, b)
            end
        end

        prev_indices = current_indices
        offset += nvertices(c)
    end

    log_info(types, "Type of points in connect_tube: ", typeof(points))

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