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
function create_circle_in_local_frame!(points::AbstractArray{T}, normals::AbstractArray{T}, center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, radius) where T

    for (i, α) in enumerate(range(0, 2*π, length=resolution+1)[1:end-1]) #alloc
        @. points[:, i] = center + radius*cos(α)*local_y + radius*sin(α)*local_z
        @. normals[:, i] = radius*cos(α)*local_y + radius*sin(α)*local_z
        normals[:, i] ./= norm(@view normals[:, i])
    end
end

"""
Returns points that form an ellipse/ a scaled circle in the plane given by local_y and local_z.
half_width is the scale factor applied to local_y, half_height is applied to local_z.  
Resolution determines the number of points. 
local_y and local_z are expected to have a length of 1. 

See also [`create_circle_in_local_frame`](@ref), [`create_rectangle_in_local_frame`](@ref).
"""
function create_ellipse_in_local_frame!(points::AbstractArray{T}, normals::AbstractArray{T}, center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, half_width, half_height) where T

    for (i, α) in enumerate(range(0, 2*π, length=resolution+1)[1:end-1]) #alloc
        @. points[:, i] = center + half_width*cos(α)*local_y + half_height*sin(α)*local_z
        @. normals[:, i] = half_width*cos(α)*local_y + half_height*sin(α)*local_z
        normals[:, i] ./= norm(@view normals[:, i])
    end
end

"""
Returns points that form a rectangle in the plane given by local_y and local_z. 
Resolution determines the number of points. 
local_y and local_z are expected to have a length of 1. 

See also [`create_circle_in_local_frame`](@ref), [`create_ellipse_in_local_frame`](@ref).
"""
function create_rectangle_in_local_frame!(points::AbstractArray{T}, normals::AbstractArray{T}, center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, half_width, half_height) where T
    if(resolution<4)
        create_ellipse_in_local_frame(points, normals, center, local_y, local_z, resolution, width, height)
    end

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
    
    half_width_vector = half_width*local_y
    half_height_vector = half_height*local_z
    a = 1
    for y in range(0, 1, A+2)[2:end-1]
        @. points[:, a] = center + half_width_vector + y*half_height_vector
        normals[:, a] .= local_y
        a+=1
    end
    @. points[:, a] = center + half_width_vector + half_height_vector
    normals[:, a] .= local_y .+ local_z
    normals[:, a] ./= norm(@view normals[:, a])
    a+=1
    for x in range(1, -1, B+2)[2:end-1]
        @. points[:, a] = center + x*half_width_vector + half_height_vector
        normals[:, a] .= local_z
        a+=1
    end
    @. points[:, a] = center - half_width_vector + half_height_vector
    normals[:, a] .= -local_y .+ local_z
    normals[:, a] ./= norm(@view normals[:, a])
    a+=1
    for y in range(1, -1, C+2)[2:end-1]
        @. points[:, a] = center - half_width_vector + y*half_height_vector
        normals[:, a] .= -local_y
        a+=1
    end
    @. points[:, a] = center - half_width_vector - half_height_vector
    normals[:, a] .= -local_y .- local_z
    normals[:, a] ./= norm(@view normals[:, a])
    a+=1
    for x in range(-1, 1, D+2)[2:end-1]
        @. points[:, a] = center + x*half_width_vector - half_height_vector
        normals[:, a] .= -local_z
        a+=1
    end
    @. points[:, a] = center + half_width_vector - half_height_vector
    normals[:, a] .= local_y .- local_z
    normals[:, a] ./= norm(@view normals[:, a])
    a+=1
    for y in range(-1, 0, E+2)[2:end-1]
        @. points[:, a] = center + half_width_vector + y*half_height_vector
        normals[:, a] .= local_y
        a+=1
    end
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
    normals = Matrix{T}(undef, 3, num_points)
    a = 1
    for m in meshes
        points[:, a:a+size(m.vertices, 2)-1] .= m.vertices
        normals[:, a:a+size(m.normals, 2)-1] .= m.normals
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

    return PlainMesh{T}(points, normals, connects, colors)
end

"""
Adds faces between circles to create the surface of a tube. 
The last two vertices are assumed to be the start- and end-cap respectively. 
"""
function add_faces_to_tube_mesh!(tube_mesh::PlainMesh{T}, resolution::Int, ncircles::Int) where {T}
    offset = 0
    connection_i = 1
    prev_indices = nothing
    # @Threads.threads 
    for index = 1:ncircles

        # add connections between circles
        prev_indices = (offset+1):(offset+resolution)
        current_indices = (offset+1):(offset+resolution)
        if(prev_indices !== nothing)
            @assert length(current_indices)==length(prev_indices)

            tube_mesh.connections[1, connection_i:connection_i+resolution-1] = current_indices'
            tube_mesh.connections[2, connection_i:connection_i+resolution-1] = prev_indices'
            tube_mesh.connections[3, connection_i:connection_i+resolution-1] = circshift(prev_indices, 1)'
            connection_i += resolution

            tube_mesh.connections[1, connection_i:connection_i+resolution-1] = current_indices'
            tube_mesh.connections[2, connection_i:connection_i+resolution-1] = circshift(prev_indices, 1)'
            tube_mesh.connections[3, connection_i:connection_i+resolution-1] = circshift(current_indices, 1)'
            connection_i += resolution
        end

        prev_indices = current_indices
        offset += resolution
    end

    start_point_index = size(tube_mesh.vertices, 2)-1
    tube_mesh.connections[1, connection_i:connection_i+resolution-1] = collect(1:resolution)'
    tube_mesh.connections[2, connection_i:connection_i+resolution-1] = circshift(1:resolution, 1)'
    tube_mesh.connections[3, connection_i:connection_i+resolution-1] = repeat([start_point_index], resolution)'
    connection_i += resolution

    tube_mesh.connections[1, connection_i:connection_i+resolution-1] = collect(start_point_index-resolution:start_point_index-1)'
    tube_mesh.connections[2, connection_i:connection_i+resolution-1] = circshift(start_point_index-resolution:start_point_index-1, 1)'
    tube_mesh.connections[3, connection_i:connection_i+resolution-1] = repeat([start_point_index+1], resolution)'
    connection_i += resolution

end