export create_circle_in_local_frame, local_frame_mesh

"Sets the color of a PlainMesh or PlainNonStdMesh."
function color!(mesh, new_color::NTuple{3, Int})
    for i=1:length(mesh.colors)
        mesh.colors[i] = new_color
    end
end

"""
Writes vertices that form a circle (in the plane given by local_y and local_z) into points (and corresponding normals into normals). 
Resolution determines the number of points. 
local_y and local_z are expected to have a length of 1. 

See also [`create_ellipse_in_local_frame`](@ref), [`create_rectangle_in_local_frame`](@ref).
"""
function create_circle_in_local_frame!(points::AbstractArray{T}, normals::AbstractArray{T}, center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, radius) where T

    for (i, α) in enumerate(range(0, 2*π, length=resolution+1)[1:end-1]) #alloc
        @. points[:, i] = center + radius*cos(α)*local_y + radius*sin(α)*local_z
        @. normals[:, i] = radius*cos(α)*local_y + radius*sin(α)*local_z
        normalize_col!(normals, i)
    end
end

"""
Writes vertices that form an ellipse/ a scaled circle (in the plane given by local_y and local_z) into points (and corresponding normals into normals). 
half_width is the scale factor applied to local_y, half_height is applied to local_z.  
Resolution determines the number of points. 
local_y and local_z are expected to have a length of 1. 

See also [`create_circle_in_local_frame`](@ref), [`create_rectangle_in_local_frame`](@ref).
"""
function create_ellipse_in_local_frame!(points::AbstractArray{T}, normals::AbstractArray{T}, center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, half_width, half_height) where T

    for (i, α) in enumerate(range(0, 2*π, length=resolution+1)[1:end-1]) #alloc
        @. points[:, i] = center + half_width*cos(α)*local_y + half_height*sin(α)*local_z
        @. normals[:, i] = half_width*cos(α)*local_y + half_height*sin(α)*local_z
        normalize_col!(normals, i)
    end
end

"""
Writes vertices that form a rectangle (in the plane given by local_y and local_z) into points (and corresponding normals into normals). 
Resolution determines the number of points. 
local_y and local_z are expected to have a length of 1. 

See also [`create_circle_in_local_frame`](@ref), [`create_ellipse_in_local_frame`](@ref).
"""
function create_rectangle_in_local_frame!(points::AbstractArray{T}, normals::AbstractArray{T}, center::AbstractArray{T}, local_y::AbstractArray{T}, local_z::AbstractArray{T}, resolution::Int, half_width, half_height) where T
    if(resolution<4)
        create_ellipse_in_local_frame!(points, normals, center, local_y, local_z, resolution, half_width, half_height)
        return
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
    normalize_col!(normals, a)
    a+=1
    for x in range(1, -1, B+2)[2:end-1]
        @. points[:, a] = center + x*half_width_vector + half_height_vector
        normals[:, a] .= local_z
        a+=1
    end
    @. points[:, a] = center - half_width_vector + half_height_vector
    normals[:, a] .= -local_y .+ local_z
    normalize_col!(normals, a)
    a+=1
    for y in range(1, -1, C+2)[2:end-1]
        @. points[:, a] = center - half_width_vector + y*half_height_vector
        normals[:, a] .= -local_y
        a+=1
    end
    @. points[:, a] = center - half_width_vector - half_height_vector
    normals[:, a] .= -local_y .- local_z
    normalize_col!(normals, a)
    a+=1
    for x in range(-1, 1, D+2)[2:end-1]
        @. points[:, a] = center + x*half_width_vector - half_height_vector
        normals[:, a] .= -local_z
        a+=1
    end
    @. points[:, a] = center + half_width_vector - half_height_vector
    normals[:, a] .= local_y .- local_z
    normalize_col!(normals, a)
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
function merge_meshes(meshes::AbstractVector{PlainMesh{T}}) where T

    points, normals, connects, colors = merge_meshes(
        map(m -> m.vertices, meshes), 
        map(m -> m.normals, meshes), 
        map(m -> m.connections, meshes), 
        map(m -> m.colors, meshes))

    return PlainMesh{T}(points, normals, connects, colors)
end

function merge_meshes(vertex_list::AbstractVector{Matrix{T}}, 
    normal_list::AbstractVector{Matrix{T}}, 
    connection_list::AbstractVector{Matrix{Int}}, 
    color_list::AbstractVector{X}) where {T, Y <: Union{String, NTuple{3, Int}}, X <: AbstractVector{Y}}

    num_points = sum(map(m -> size(m, 2), vertex_list))
    num_connects = sum(map(m -> size(m, 2), connection_list))
    
    points = Matrix{T}(undef, 3, num_points)
    normals = Matrix{T}(undef, 3, num_points)
    connects = Matrix{Int}(undef, 3, num_connects)
    colors = Vector{Y}(undef, num_points)

    point_count = 0
    connection_count = 0

    for (v, n, con, col) in zip(vertex_list, normal_list, connection_list, color_list)
        num = size(v, 2)
        points[:, point_count+1:point_count+num] .= v
        normals[:, point_count+1:point_count+num] .= n
        colors[point_count+1:point_count+num] .= col

        connects[:, connection_count+1:connection_count+size(con, 2)] .= con .+ point_count

        point_count += num
        connection_count += size(con, 2)
    end

    return points, normals, connects, colors
end

"""
Adds faces between circles to create the surface of a tube. 
The last two vertices are assumed to be the start- and end-cap respectively. 
The resolution is the number of vertices per circle/cross_section. 
"""
function add_faces_to_tube_mesh!(tube_mesh::PlainMesh{T}, resolution::Int, ncircles::Int) where {T}

    if ncircles<=0 || size(tube_mesh.vertices, 2)<=0
        return
    end

    num_faces = (ncircles-1)*resolution*2 + 2*resolution # connections between neighbor circles + connections to the end points
    tube_mesh.connections = Matrix{Int}(undef, 3, num_faces)

    offset = 0
    connection_i = 1
    prev_indices = nothing
    # @Threads.threads 
    for index = 1:ncircles

        # add connections between circles
        current_indices = (offset+1):(offset+resolution)
        shift_buffer = Vector{Int}(undef, resolution)
        if(prev_indices !== nothing)
            @assert length(current_indices)==length(prev_indices)

            tube_mesh.connections[1, connection_i:connection_i+resolution-1] = current_indices'
            tube_mesh.connections[2, connection_i:connection_i+resolution-1] = prev_indices'
            circshift!(shift_buffer, prev_indices, 1)
            tube_mesh.connections[3, connection_i:connection_i+resolution-1] = shift_buffer'
            connection_i += resolution

            tube_mesh.connections[1, connection_i:connection_i+resolution-1] = current_indices'
            tube_mesh.connections[2, connection_i:connection_i+resolution-1] = shift_buffer'
            circshift!(shift_buffer, current_indices, 1)
            tube_mesh.connections[3, connection_i:connection_i+resolution-1] = shift_buffer'
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