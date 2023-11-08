export create_circle_in_local_frame, local_frame_mesh

function color!(mesh, new_color)
    mesh.colors = repeat([new_color], size(mesh.vertices, 2))
end

function translate!(mesh, vec::AbstractArray{T}) where T <: AbstractFloat 
    mesh.vertices .+= vec
end

# Assumes the object in the origin and rotates it, so that its new local z-axis points to direction. 
function rotate_in_direction!(mesh::Union{PlainMesh, PlainNonStdMesh}, direction::AbstractArray{T}) where T<:AbstractFloat

    if(direction[1]!=0)
        rotationAxis = Vec(-direction[2]/direction[1], T(1), T(0)) * -sign(direction[1])
    elseif(direction[2]!=0)
        rotationAxis = Vec(T(1), -direction[1]/direction[2], T(0)) * sign(direction[2]) # rotation axis is perpendicular to direction and lies in the xy-plane
    else
        return
    end
    rotationAngle = Meshes.∠(Vec{3,T}(direction...), Vec(T(0), T(0), T(1)))
    mesh.vertices = (mesh.vertices' * AngleAxis(rotationAngle, rotationAxis...))'
end

function create_circle_in_local_frame(center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, radius) where T<:AbstractFloat
    points = Matrix{T}(undef, 3, resolution)

    for (i, α) in enumerate(collect(range(0, 2*π, length=resolution+1))[1:end-1])
        points[:, i] = center + radius*cos(α)*local_y + radius*sin(α)*local_z
    end
    return points
end

function create_ellipse_in_local_frame(center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, half_width, half_height) where T<:AbstractFloat
    points = Matrix{T}(undef, 3, resolution)

    for (i, α) in enumerate(collect(range(0, 2*π, length=resolution+1))[1:end-1])
        points[:, i] = center + half_width*cos(α)*local_y + half_height*sin(α)*local_z
    end
    return points
end

function create_rectangle_in_local_frame(center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, half_width, half_height) where T<:AbstractFloat
    if(resolution<4)
        return create_ellipse_in_local_frame(center, local_y, local_z, resolution, width, height)
    end

    points = Matrix{T}(undef, 3, resolution)

    # each corner has a point
    # the remaining points are split into sections A to E
    remaining_points = resolution-4
    ratio = 1/6
    short_sides = Int(ceil(remaining_points*ratio))
    long_sides = remaining_points - short_sides
    
    B = Int(floor(long_sides/2))
    D = long_sides-B
    C = Int(floor(short_sides/2))
    AE = short_sides-C
    A = Int(floor(AE/2))
    E = AE-A
    
    half_width = half_width*local_y
    half_height = half_height*local_z
    a = 1
    for y in collect(range(0, 1, A+2))[2:end-1]
        points[:, a] = center + half_width + y*half_height
        a+=1
    end
    points[:, a] = center + half_width + half_height
    a+=1
    for x in collect(range(1, -1, B+2))[2:end-1]
        points[:, a] = center + x*half_width + half_height
        a+=1
    end
    points[:, a] = center - half_width + half_height
    a+=1
    for y in collect(range(1, -1, C+2))[2:end-1]
        points[:, a] = center - half_width + y*half_height
        a+=1
    end
    points[:, a] = center - half_width - half_height
    a+=1
    for x in collect(range(-1, 1, D+2))[2:end-1]
        points[:, a] = center + x*half_width - half_height
        a+=1
    end
    points[:, a] = center + half_width - half_height
    a+=1
    for y in collect(range(-1, 0, E+2))[2:end-1]
        points[:, a] = center + half_width + y*half_height
        a+=1
    end
    return points
end

function local_frame_mesh(local_zero, local_x, local_y, local_z)
    x = ColoredMesh(discretize(CylinderSurface(Tuple(local_zero), Tuple(local_zero + local_x), 0.01), RegularDiscretization(6)), (255, 0, 0))
    y = ColoredMesh(discretize(CylinderSurface(Tuple(local_zero), Tuple(local_zero + local_y), 0.01), RegularDiscretization(6)), (0, 255, 0))
    z = ColoredMesh(discretize(CylinderSurface(Tuple(local_zero), Tuple(local_zero + local_z), 0.01), RegularDiscretization(6)), (0, 0, 255))

    a = merge(x, y)
    return merge(a, z)
end

function local_arrow_mesh(center, vector, color)
    return ColoredMesh(discretize(CylinderSurface(0.01, Segment([Tuple(center), Tuple(center + vector)])), RegularDiscretization(6)), color)
end



        

# Merges multiple mesh objects into one. The geometry (vertices/edges/faces) does not change. 
function merge_multiple_meshes(meshes::AbstractVector{PlainMesh{T}}) where {T}
    num_points = sum(map(m -> size(m.vertices, 2), meshes))
    num_connects = sum(map(m -> size(m.connections, 2), meshes))

    points = hcat(map(m->m.vertices, meshes)...)

    colors = Array{NTuple{3, Int}}(undef, num_points)
    connects = Array{Int, 2}(undef, 3, num_connects)

    point_offset = 0
    connect_offset = 0
    for m in meshes
        point_len = size(m.vertices, 2)
        connect_len = size(m.connections, 2)

        
        colors[point_offset+1 : point_offset+point_len] = m.colors
        

        connects[:, connect_offset+1 : connect_offset+connect_len] = m.connections .+ point_offset

        point_offset += point_len
        connect_offset += connect_len
    end

    return PlainMesh(points, connects, colors)
end

# Adds faces between circles to create the surface of a tube. 
function connect_circles_to_tube(circles::AbstractVector{PlainNonStdMesh{T}}, endpoints = nothing) where {T}

    shiftSum = 0

    # collect all points
    if(endpoints===nothing)
        points = hcat(map(m->m.vertices, circles)...)
        colors = vcat(map(c -> c.colors, circles)...)
    else
        @assert length(endpoints)==2
        points = hcat(map(m->m.vertices, circles)..., endpoints...)
        colors = vcat(map(c -> c.colors, circles)..., circles[1].colors[1], circles[end].colors[1])
    end

    resolution = size(circles[1].vertices, 2)
    connections = Array{Int, 2}(undef, 3, (length(circles)-1)*(2*resolution) + (endpoints===nothing ? 0 : 2*resolution))
    
    offset = 0
    connection_i = 1
    prev_indices = nothing
    color_count = 0
    for c in circles

        # add connections between circles
        current_indices = (offset+1):(offset+nvertices(c))
        if(prev_indices !== nothing)
            @assert length(current_indices)==length(prev_indices)

            connections[1, connection_i:connection_i+resolution-1] = current_indices'
            connections[2, connection_i:connection_i+resolution-1] = prev_indices'
            connections[3, connection_i:connection_i+resolution-1] = circshift(prev_indices, 1)'
            connection_i += resolution

            connections[1, connection_i:connection_i+resolution-1] = current_indices'
            connections[2, connection_i:connection_i+resolution-1] = circshift(prev_indices, 1)'
            connections[3, connection_i:connection_i+resolution-1] = circshift(current_indices, 1)'
            connection_i += resolution
        end

        prev_indices = current_indices
        offset += nvertices(c)
    end

    if(endpoints!==nothing) # create ends
        start_point_index = size(points, 2)-1
        connections[1, connection_i:connection_i+resolution-1] = collect(1:resolution)'
        connections[2, connection_i:connection_i+resolution-1] = circshift(1:resolution, 1)'
        connections[3, connection_i:connection_i+resolution-1] = repeat([start_point_index], resolution)'
        connection_i += resolution

        connections[1, connection_i:connection_i+resolution-1] = collect(start_point_index-resolution:start_point_index-1)'
        connections[2, connection_i:connection_i+resolution-1] = circshift(start_point_index-resolution:start_point_index-1, 1)'
        connections[3, connection_i:connection_i+resolution-1] = repeat([start_point_index+1], resolution)'
        connection_i += resolution
    end


    log_info(types, "Type of points in connect_tube: ", typeof(points))

    log_info(frame_rotation, "Total needed rotation: ", shiftSum)

    
    return PlainMesh(points, connections, colors)
end

function determine_offset(p1, p2, circle_points)

    distances1 = map(cp -> norm(p1 .- cp), eachcol(circle_points))
    nearest1 = argmin(distances1)

    distances2 = map(cp -> norm(p2 .- cp), eachcol(circle_points))
    distances2[nearest1] = max(distances2...)+1
    nearest2 = argmin(distances2)

    a = mod(nearest2-nearest1, size(circle_points, 2))
  
    if(a>size(circle_points, 2)/2)

        log_info(circle_index_correction, nearest1, " ", nearest2)
        return nearest1-1, true
    else 
        return nearest1-1, false
    end
end