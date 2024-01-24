import Base: push!, iterate, IteratorSize, eltype, isdone
"""
A HashGrid contains 3D points and can be used to access points in a neighborhood efficiently. 
The 3D points can be inside of objects of type U. 
Type T determines the precision of the grid(float32/Flaot64). 
"""
mutable struct HashGrid{T, U}
    data::Vector{U}
    grid::Array{Vector{Int}, 3}
    position_access_function
    origin::Vector{T}
    box_size::Vector{T}

    function HashGrid{T, U}(bounding_box::Matrix{T}, box_size::Vector{T}, position_access_function, estimated_n_entries=0) where {T, U}
        if size(bounding_box)!=(3, 2)
            throw(ArgumentError("The bounding box should be a 3x2 matrix with the minimal values in the first and the maximal values in the second column. "))
        end
        if length(box_size)!=3
            throw(ArgumentError("The box size should be a vector with 3 elements. "))
        end
        for i=1:3
            if box_size[i]<=0
                throw(ArgumentError("The box size has to be positive in all dimensions. "))
            end
            if !(bounding_box[:, 1]<=bounding_box[:, 2])
                throw(ArgumentError("The max value in the bounding box needs to be greater or equal than the min value. "))
            end
        end
        estimated_n_entries = max(estimated_n_entries, 0)

        bounding_box[:, 2] .+= 0.1 # the real bounding box needs to be a bit larger to include the upper end of the bounding box

        n_boxes = @. convert(Int, (ceil( (bounding_box[:, 2]-bounding_box[:, 1])/box_size )))
        data = Vector{U}()
        sizehint!(data, estimated_n_entries)
        grid = Array{Vector{Int}}(undef, n_boxes...)
        for i=eachindex(grid)
            grid[i] = Vector{Int}()
        end
        new{T, U}(data, grid, position_access_function, bounding_box[:, 1], box_size)
    end
end

struct HashGridIterationHelper{T, U}
    grid::HashGrid{T, U}
    world_coordinates::Vector{T}
    grid_index::Vector{Int}
    max_distance::T

    function HashGridIterationHelper{T, U}(grid::HashGrid{T, U}, world_coordinates::Vector{T}, max_distance::T) where {T, U}
        index = world_to_index(grid, world_coordinates)
        return new{T, U}(grid, world_coordinates, index, max_distance)
    end
end



function world_to_index(grid, world_coordinates)
    index_list = [0, 0, 0]
    for i=1:3
        index_list[i] = convert(Int, floor((world_coordinates[i]-grid.origin[i])/grid.box_size[i])) + 1
        if index_list[i] < 1 || index_list[i] > size(grid.grid, i)
            throw(DomainError(world_coordinates, "Coordinates were outside of the grid. "))
        end
    end
    return index_list
end

function Base.push!(grid::HashGrid{T, U}, entries::Vararg{U}) where {T, U}
    for entry in entries
        push!(grid.data, entry)
        index_list = world_to_index(grid, grid.position_access_function(entry))
        push!(grid.grid[index_list...], length(grid.data))
    end
end

function iteration_recursion(grid, central_index, delta_list, layer)
    if layer<1
        return nothing
    end
    
    if delta_list[layer]==1 || central_index[layer]+delta_list[layer] == size(grid.grid, layer)
        delta_list[layer] = central_index[layer]==1 ? 0 : -1
        return iteration_recursion(grid, central_index, delta_list, layer-1)
    else
        delta_list[layer] += 1
        return delta_list
    end
end

function iterate(item::HashGridIterationHelper)
    # the state is a tuple (delta_i, delta_j, delta_k, local_index)
    # the delta values are in [-1, 1] and select the current neighbor box in the hashgrid. They are usually initialized with -1, but when item.grid_index is at the left border of the grid, it is 0. 
    # local_index is either -1 or an index for accessing the index-vector in the current box

    return Base.iterate(item, (item.grid_index[1]==1 ? 0 : -1, item.grid_index[2]==1 ? 0 : -1, (item.grid_index[3]==1 ? 0 : -1)-1, -1))
end

function iterate(item::HashGridIterationHelper, (i, j, k, local_next_index))
    if(local_next_index==-1)
        log_info(hash_grid, "no local index")
        new_ijk = iteration_recursion(item.grid, item.grid_index, [i, j, k], 3)
        if(new_ijk===nothing)
            log_info(hash_grid, "no data left :(")
            return nothing
        else
            log_info(hash_grid, "new box", new_ijk)
            if length(item.grid.grid[(item.grid_index .+ new_ijk)...]) > 0
                return Base.iterate(item, (new_ijk..., 1))
            else
                return Base.iterate(item, (new_ijk..., -1))
            end
        end

    else
        data_index = item.grid.grid[(item.grid_index .+ (i, j, k))...][local_next_index]
        value = item.grid.data[data_index]
        local_next_index += 1
        if local_next_index>length(item.grid.grid[(item.grid_index .+ (i, j, k))...])
            local_next_index = -1
        end
        dist = norm(item.world_coordinates .- item.grid.position_access_function(value))
        if(dist>item.max_distance) # skip entries that are in adjacent boxes but still too far away
            log_info(hash_grid, "too far away ($dist)", item.world_coordinates, item.grid.position_access_function(value))
            return Base.iterate(item, (i, j, k, local_next_index))
        else
            log_info(hash_grid, "found point", value)
            return value, (i, j, k, local_next_index)
        end
    end
end

Base.IteratorSize(::Type{HashGridIterationHelper{T, U}}) where {T, U} = Base.SizeUnknown()
Base.eltype(::Type{HashGridIterationHelper{T, U}}) where {T, U} = U
function Base.isdone(item::HashGridIterationHelper, (i, j, k, local_next_index))
    if local_next_index==-1
        while true
            new_ijk = iteration_recursion(item.grid, item.grid_index, [i, j, k], 3)
            if(new_ijk===nothing)
                return true
            else
                if length(item.grid.grid[(item.grid_index .+ new_ijk)...]) > 0
                    return false
                else
                    i = new_ijk[1]
                    j = new_ijk[2]
                    k = new_ijk[3]
                end
            end
        end
    else
        return false
    end
end


function each_neighbor(grid::HashGrid{T, U}, world_coordinates::AbstractVector{T}, max_distance::T) where {T, U}# TODO exclude self
    return HashGridIterationHelper{T,U}(grid, world_coordinates, max_distance)
end




# get_indices|entries_near_indexed_data|world_pos()