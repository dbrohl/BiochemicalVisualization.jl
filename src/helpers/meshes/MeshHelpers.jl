export create_circle_in_local_frame, local_frame_mesh

"Sets the color of a PlainMesh or PlainNonStdMesh."
function color!(mesh, new_color::NTuple{3, Int})
    for i=1:length(mesh.colors)
        mesh.colors[i] = new_color
    end
end

"""
Returns points that form a circle in the plane given by local_y and local_z. 
Resolution determines the number of points. 
local_y and local_z are expected to have a length of 1. 

See also [`create_ellipse_in_local_frame`](@ref), [`create_rectangle_in_local_frame`](@ref).
"""
function create_circle_in_local_frame(center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, radius) where T
    points = Matrix{T}(undef, 3, resolution)

    for (i, α) in enumerate(range(0, 2*π, length=resolution+1)[1:end-1]) #alloc
        @. points[:, i] = center + radius*cos(α)*local_y + radius*sin(α)*local_z
    end
    return points
end

"""
Returns points that form an ellipse/ a scaled circle in the plane given by local_y and local_z.
half_width is the scale factor applied to local_y, half_height is applied to local_z.  
Resolution determines the number of points. 
local_y and local_z are expected to have a length of 1. 

See also [`create_circle_in_local_frame`](@ref), [`create_rectangle_in_local_frame`](@ref).
"""
function create_ellipse_in_local_frame(center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, half_width, half_height) where T
    points = Matrix{T}(undef, 3, resolution)

    for (i, α) in enumerate(range(0, 2*π, length=resolution+1)[1:end-1]) #alloc
        @. points[:, i] = center + half_width*cos(α)*local_y + half_height*sin(α)*local_z
    end
    return points
end

"""
Returns points that form a rectangle in the plane given by local_y and local_z. 
Resolution determines the number of points. 
local_y and local_z are expected to have a length of 1. 

See also [`create_circle_in_local_frame`](@ref), [`create_ellipse_in_local_frame`](@ref).
"""
function create_rectangle_in_local_frame(center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, half_width, half_height) where T
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
    for y in range(0, 1, A+2)[2:end-1]
        @. points[:, a] = center + half_width + y*half_height
        a+=1
    end
    @. points[:, a] = center + half_width + half_height
    a+=1
    for x in range(1, -1, B+2)[2:end-1]
        @. points[:, a] = center + x*half_width + half_height
        a+=1
    end
    @. points[:, a] = center - half_width + half_height
    a+=1
    for y in range(1, -1, C+2)[2:end-1]
        @. points[:, a] = center - half_width + y*half_height
        a+=1
    end
    @. points[:, a] = center - half_width - half_height
    a+=1
    for x in range(-1, 1, D+2)[2:end-1]
        @. points[:, a] = center + x*half_width - half_height
        a+=1
    end
    @. points[:, a] = center + half_width - half_height
    a+=1
    for y in range(-1, 0, E+2)[2:end-1]
        @. points[:, a] = center + half_width + y*half_height
        a+=1
    end
    return points
end

"Creates a ColoredMesh representing a red, green and blue coordinate system. "
function local_frame_mesh(local_zero, local_x, local_y, local_z)
    x = ColoredMesh(Meshes.discretize(Meshes.CylinderSurface(Tuple(local_zero), Tuple(local_zero + local_x), 0.01), Meshes.RegularDiscretization(6)), (255, 0, 0))
    y = ColoredMesh(Meshes.discretize(Meshes.CylinderSurface(Tuple(local_zero), Tuple(local_zero + local_y), 0.01), Meshes.RegularDiscretization(6)), (0, 255, 0))
    z = ColoredMesh(Meshes.discretize(Meshes.CylinderSurface(Tuple(local_zero), Tuple(local_zero + local_z), 0.01), Meshes.RegularDiscretization(6)), (0, 0, 255))

    a = merge(x, y)
    return merge(a, z)
end

"Creates a ColoredMesh that points from the center into the direction of vector. "
function local_arrow_mesh(center, vector, color)
    return ColoredMesh(Meshes.discretize(Meshes.CylinderSurface(Tuple(center), Tuple(center + vector), 0.01), Meshes.RegularDiscretization(6)), color)
end

"Merges multiple mesh objects into one. The geometry (vertices/edges/faces) does not change. "
function merge_multiple_meshes(meshes::AbstractVector{PlainMesh{T}}) where T
    num_points = sum(map(m -> size(m.vertices, 2), meshes))
    num_connects = sum(map(m -> size(m.connections, 2), meshes))

    points = Matrix{T}(undef, 3, num_points)
    a = 1
    for m in meshes
        points[:, a:a+size(m.vertices, 2)-1] .= m.vertices
        a += size(m.vertices, 2)
    end

    colors = Array{NTuple{3, Int}}(undef, num_points)
    connects = Array{Int, 2}(undef, 3, num_connects)

    point_offset = 0
    connect_offset = 0
    for m in meshes
        point_len = size(m.vertices, 2)
        connect_len = size(m.connections, 2)

        
        colors[point_offset+1 : point_offset+point_len] = m.colors
        

        connects[:, connect_offset+1 : connect_offset+connect_len] .= m.connections .+ point_offset

        point_offset += point_len
        connect_offset += connect_len
    end

    return PlainMesh{T}(points, connects, colors)
end

"""
Adds faces between circles to create the surface of a tube. 
When endpoints is passed to the function, the tube will be closed (by a flat plane) at the start and the end. 
"""
function connect_circles_to_tube(circles::AbstractVector{PlainNonStdMesh{T}}, endpoints::Union{Nothing, NTuple{2, AbstractVector{T}}} = nothing) where {T}

    if(length(circles)==0)
        return PlainMesh(Matrix{Int}(undef, 3, 0), Matrix{Int}(undef, 3, 0), Vector{NTuple{3, Int}}())
    end

    resolution = size(circles[1].vertices, 2)

    # collect all points
    matrix_length = (length(circles) * resolution) + (endpoints===nothing ? 0 : 2)
    points = Matrix{T}(undef, 3, matrix_length)
    colors = Vector{NTuple{3, Int}}(undef, matrix_length)
    for (i, circle) in enumerate(circles)
        points[:, (i-1)*resolution+1:i*resolution] = circle.vertices
        colors[(i-1)*resolution+1:i*resolution] = circle.colors
    end
    if(endpoints!==nothing)
        points[:, end-1] = endpoints[1]
        colors[end-1] = circles[1].colors[1]
        points[:, end] = endpoints[2]
        colors[end] = circles[end].colors[1]
    end

    
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

    
    return PlainMesh{T}(points, connections, colors)
end