export prepare_backbone_model, generate_chain_mesh, generate_geometry_at_point, benchmark_method # TODO remove

function benchmark_method(ac::System{T}, config::BackboneConfig) where T
    mesh = prepare_backbone_model(ac, config)
    representation = Representation(mesh) 
    return representation
end

"""
Inserts elem into array. Assuming that array is sorted ascending, it will be still sorted after the insertion. 
"""
function insert_sorted!(array, elem)
    index = searchsortedfirst(array, elem)
    insert!(array, index, elem)
end

"""
Loops over a the list sample_indices and counts observed secondary structures. 
- sample_indices contains keys for the residue_info_dict. 
- residue_info_dict contains tuple-values where the second element is the secondary structure
"""
function get_ss_count(sample_indices, residue_info_dict)
    ss_count = Dict{BiochemicalAlgorithms.SecondaryStructure.T, Int}(
        BiochemicalAlgorithms.SecondaryStructure.NONE => 0, 
        BiochemicalAlgorithms.SecondaryStructure.HELIX => 0, 
        BiochemicalAlgorithms.SecondaryStructure.SHEET => 0
    )

    prev_ss = nothing
    for index in sample_indices
        curr_ss = residue_info_dict[index][2]
        if(prev_ss!=curr_ss)
            ss_count[curr_ss]+=1
            prev_ss = curr_ss
        end
    end
    return ss_count
end

"""
Finds amino acids that are the end of a beta sheet and should be displayed as an arrowhead.
"""
function compute_frame_widths(fragment_list::Vector{Fragment{T}}, sample_to_residue_indices, residue_info_dict::Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}}) where T
    rectangle_widths::Vector{T} = []
    arrow_starts::Vector{Int} = []
    arrow_frame_indices::Vector{Int} = []

    filtered_fragment_list = filter(f -> is_amino_acid(f), fragment_list)
    # determine orientation of looping (arrows point towards carboxyl-end)
    if((Symbol(:N_TERMINAL) ∈ filtered_fragment_list[1].flags) 
        && (Symbol(:C_TERMINAL) ∈ filtered_fragment_list[end].flags))
        n_to_c = true
    elseif((Symbol(:N_TERMINAL) ∈ filtered_fragment_list[end].flags) 
        && (Symbol(:C_TERMINAL) ∈ filtered_fragment_list[1].flags))
        n_to_c = false
    else
        throw(ErrorException("c and n terminal are not included in flags (chain $chain_num)"))
    end
    log_info(misc, "n_to_c: $n_to_c")
    
    # for each sheet: find the end and save the index
    prev_idx = nothing
    prev_ss = nothing
    arrow_fragment_indices = []


    fragment_idx_list_from_frames = unique(sample_to_residue_indices)
    for current_fragment_idx in fragment_idx_list_from_frames
        current_ss = residue_info_dict[current_fragment_idx][2]

        if(n_to_c && prev_ss!==nothing
            && prev_ss==BiochemicalAlgorithms.SecondaryStructure.SHEET 
            && current_ss!=BiochemicalAlgorithms.SecondaryStructure.SHEET)
            push!(arrow_fragment_indices, prev_idx)
        end

        if(!n_to_c 
            && prev_ss!=BiochemicalAlgorithms.SecondaryStructure.SHEET 
            && current_ss==BiochemicalAlgorithms.SecondaryStructure.SHEET)
            push!(arrow_fragment_indices, current_fragment_idx)
        end
        prev_idx = current_fragment_idx
        prev_ss = current_ss
    end


    if(n_to_c && prev_ss==BiochemicalAlgorithms.SecondaryStructure.SHEET) # end last sheet
        push!(arrow_fragment_indices, prev_idx)
    end

    # find corresponding frames and store the resulting width of the rectangle at that frame
    frames_in_residue_count = 0
    prev_res_idx = sample_to_residue_indices[1]

    a = 1
    while(a<=length(sample_to_residue_indices))
        
        res_idx = sample_to_residue_indices[a]
        is_last_frame = a==length(sample_to_residue_indices)
        if(res_idx!=prev_res_idx || is_last_frame) # end of the previous residue
            if(is_last_frame)
                frames_in_residue_count += 1
            end
            if(prev_res_idx ∈ arrow_fragment_indices || (is_last_frame && res_idx ∈ arrow_fragment_indices))
                num_uniform = Int(round(frames_in_residue_count/3))
                num_arrow = frames_in_residue_count - num_uniform
                uniforms = repeat([1], num_uniform)
                arrow = collect(range(1.5, 0, num_arrow))
                if(n_to_c)
                    append!(arrow_frame_indices, length(rectangle_widths)+num_uniform+1:length(rectangle_widths)+num_uniform+num_arrow)
                    push!(arrow_starts, length(rectangle_widths)+num_uniform+1)
                    append!(rectangle_widths, uniforms, arrow)
                else
                    append!(arrow_frame_indices, length(rectangle_widths)+1:length(rectangle_widths)+num_arrow)
                    push!(arrow_starts, length(rectangle_widths)+num_arrow)
                    append!(rectangle_widths, reverse!(arrow), uniforms)
                end
            else
                append!(rectangle_widths, repeat([1], frames_in_residue_count))
            end

            # start new residue
            frames_in_residue_count = 0
            prev_res_idx = res_idx
        end
        frames_in_residue_count += 1
        a += 1
    end

    return rectangle_widths, arrow_starts, arrow_frame_indices, n_to_c
end

"""
Creates a single circle/ellipse/rectangle and writes the result into result_mesh. 

# Arguments
- `result_mesh::PlainMesh{T}`: A preallocated struct that will contain the result. 
- `result_mesh_index::Int`: The created frame is the result_mesh_index-th frame along the spline. This parameter determines where in result_mesh's arrays, data is inserted. 
- `point::AbstractVector{T}`: The position of the frame
- `tangent::AbstractVector{T`
- `normal::AbstractVector{T}`
- `binormal::AbstractVector{T}`
- `residue_info_dict::Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}`: Maps indices of residues to tuples (3 letter amino acid name, secondary structure)
- `linked_residue_idx::Union{Nothing,Int}`: used to look up the resiude/secondary structure in residue_info_dict
- `frame_config::Union{Nothing,Tuple{Bool, Bool, Int, Int}}`: If the point is a spline_point, frame_config is nothing. If the point is a transition_point, it has both properties from the previous and the following point. These are stored as (XXXX)
- `rectangle_width::T`: Normally 1.0, except for frames that are part of an arrow head
- `fixed_color::Union{NTuple{3, Int}, Nothing}`: contains the fixed color if there is one (e. g. fixed color for the whole chain), otherwise nothing
- `config::BackboneConfig{T}`
"""
function generate_geometry_at_point!(
    result_mesh::PlainMesh{T},
    result_mesh_index::Int, 
    point::AbstractVector{T}, 
    tangent::AbstractVector{T},
    normal::AbstractVector{T}, 
    binormal::AbstractVector{T}, 
    residue_info_dict::Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}},
    linked_residue_idx::Union{Nothing,Int}, 
    frame_config::Union{Nothing,Tuple{Bool, Bool, Int, Int}}, 
    rectangle_width::T, 
    fixed_color::Union{NTuple{3, Int}, Nothing}, 
    config::BackboneConfig{T}) where T
    # generate cross-section vertices
    start_index = (result_mesh_index-1)*config.resolution+1
    end_index = result_mesh_index*config.resolution
    if(config.backbone_type==BackboneType.BACKBONE)
        create_circle_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution, config.stick_radius)
    elseif(config.backbone_type==BackboneType.RIBBON)
        create_ellipse_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution, T(3)*config.stick_radius, config.stick_radius)
    elseif(config.backbone_type==BackboneType.CARTOON)
        shortcut = false
        if(linked_residue_idx===nothing)
            if((frame_config[1] == frame_config[2]))
                result_mesh.vertices[:, start_index:end_index] = stack(repeat([point], config.resolution)) # only 1 instead of resolution vertices would suffice, but then connect_circles_to_tube has to be modified
                result_mesh.normals[:, start_index:end_index] = stack(repeat([tangent], config.resolution))
                shortcut = true
            else
                residue_idx = frame_config[3]
            end
        else
            residue_idx = linked_residue_idx
        end

        if(!shortcut)
            structure = residue_info_dict[residue_idx][2]
            if(structure==BiochemicalAlgorithms.SecondaryStructure.NONE)
                create_circle_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution, config.stick_radius)
            elseif(structure==BiochemicalAlgorithms.SecondaryStructure.HELIX)
                create_ellipse_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution, T(3)*config.stick_radius, T(1.5)*config.stick_radius)
            elseif(structure==BiochemicalAlgorithms.SecondaryStructure.SHEET)
                create_rectangle_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution, T(3)*config.stick_radius * rectangle_width, T(0.5)*config.stick_radius)
            end
        end
    end

    # color
    color = nothing
    if(fixed_color!==nothing)
        color = fixed_color
    elseif(config.color==Color.SECONDARY_STRUCTURE)
        if(linked_residue_idx===nothing)
            residue_idx = frame_config[4]
        else
            residue_idx = linked_residue_idx
        end
        structure = residue_info_dict[residue_idx][2]
        color = SS_COLORS[structure]
    elseif(config.color==Color.RESIDUE)
        if(linked_residue_idx===nothing)
            residue_idx = frame_config[4]
        else
            residue_idx = linked_residue_idx
        end
        color = AA_COLORS[residue_info_dict[residue_idx][1]]
    end
    for i=start_index:end_index
        result_mesh.colors[i] = color
    end
end

"""
Generates a PlainMesh for chain. The config should be checked previously. 
When the whole mesh should have a uniform color, it can be passed as fixed_color. 
"""
function prepare_backbone_model(chain::Chain{T}, config::BackboneConfig{T}, fixed_color::Union{Nothing, NTuple{3, Int}} = nothing) where {T<:Real}
    if(fixed_color===nothing && (config.color==Color.UNIFORM || config.color==Color.CHAIN))
        fixed_color = (0, 0, 255)
    end

    fragment_list = fragments(chain) #TODO fragments/eachfragment

    # construct spline
    local spline
    if(config.spline==Spline.CATMULL_ROM)
        spline = CatmullRom(chain, config.control_point_strategy)
    elseif(config.spline==Spline.CUBIC_B)
        spline = CubicB(chain, config.control_point_strategy)
    end

    vertices_per_unit = T(0.4 * config.resolution / (2*π*config.stick_radius))
    # sample along spline
    spline_points, sample_to_residue_indices::Vector{Union{Int, Nothing}} = calculate_points(spline, vertices_per_unit) #alloc
    velocities = calculate_velocities(spline, vertices_per_unit)
    
    # construct local frames
    local q::Matrix{T}
    local r::Matrix{T}
    local s::Matrix{T}
    if(config.frame==Frame.RMF)
        q, r, s = rmf(spline_points, velocities)
    elseif(config.frame==Frame.SECOND_SPLINE)
        second_spline_points = calculate_minor_points(spline, vertices_per_unit)
        q, r, s = frames_from_two_splines(spline_points, velocities, second_spline_points)
    end


    # when secondary structure is displayed, additional points are necessary
    ss_count = get_ss_count(sample_to_residue_indices, spline.residue_info_dict)
    num_transition_points = 0
    if(config.backbone_type==BackboneType.CARTOON)
        num_transition_points += 2*(ss_count[BiochemicalAlgorithms.SecondaryStructure.NONE]+ss_count[BiochemicalAlgorithms.SecondaryStructure.HELIX]+ss_count[BiochemicalAlgorithms.SecondaryStructure.SHEET]-1) # changes in secondary structure
        num_transition_points += ss_count[BiochemicalAlgorithms.SecondaryStructure.SHEET] # start of arrow heads
    end

    # The additional points will not be inserted into spline_points because of efficient memory allocation. 
    # Instead we store them (with related additional data) in a separate vector
    # The tuples contain (position, q, r, s, rectangle_width, residue_index, (is_small_to_large_transition, first_transition_part, residue_for_structure, residue_for_color))
    transition_data = Vector{Tuple{Vector{T}, Vector{T}, Vector{T}, Vector{T}, T, Union{Int, Nothing}, Union{Nothing, Tuple{Bool, Bool, Int, Int}}}}(undef, num_transition_points)
    # we store before which points of spline_points each transition_point should be inserted. 
    # the array will be sorted in ascending order to "merge" it with spline_points easier
    transition_insertion_indices = Vector{Int}(undef, num_transition_points)

    fixed_indices::Vector{Int} = [] # collection of all spline_points that should not be removed by filtering 

    arrow_insert_indices = collect(1:ss_count[BiochemicalAlgorithms.SecondaryStructure.SHEET]) # stores indices that are left out for arrow transition points later
    # add frames when secondary structure changes and at arrow head starts
    if(config.backbone_type==BackboneType.CARTOON)
        # arrows part 1
        rectangle_widths, arrow_starts, arrow_point_indices, n_to_c = compute_frame_widths(fragment_list, sample_to_residue_indices, spline.residue_info_dict)
        append!(fixed_indices, arrow_point_indices) # the arrow part with changing frame widths should not be discarded by filter methods

        # ss changes
        prev_res_idx = sample_to_residue_indices[1]
        a = 1
        b = 1

        if spline.residue_info_dict[prev_res_idx][2]==BiochemicalAlgorithms.SecondaryStructure.SHEET # detect sheets at the beginning
            arrow_insert_indices[a] = b
            a+=1
            b+=1
        end

        for i=eachindex(sample_to_residue_indices)
            
            res_idx = sample_to_residue_indices[i]
            prev_ss = spline.residue_info_dict[prev_res_idx][2]
            curr_ss = spline.residue_info_dict[res_idx][2]

        
            if(res_idx!=prev_res_idx && prev_ss!=curr_ss) # ss change
                small_to_large = ((n_to_c && (prev_ss==BiochemicalAlgorithms.SecondaryStructure.NONE || prev_ss==BiochemicalAlgorithms.SecondaryStructure.SHEET)) || (!n_to_c && curr_ss==BiochemicalAlgorithms.SecondaryStructure.HELIX))
                if(small_to_large)
                    insertion_idx = i-1
                    transition_data[b] = (spline_points[:, insertion_idx], q[:, insertion_idx], r[:, insertion_idx], s[:, insertion_idx], rectangle_widths[insertion_idx], nothing, (small_to_large, true, prev_res_idx, res_idx))
                    transition_insertion_indices[b] = insertion_idx+1
                    b+=1
                    transition_data[b] = (spline_points[:, insertion_idx], q[:, insertion_idx], r[:, insertion_idx], s[:, insertion_idx], rectangle_widths[insertion_idx], nothing, (small_to_large, false, res_idx, res_idx))
                    transition_insertion_indices[b] = insertion_idx+1
                    b+=1
                    
                else
                    insertion_idx = i
                    transition_data[b] = (spline_points[:, insertion_idx], q[:, insertion_idx], r[:, insertion_idx], s[:, insertion_idx], rectangle_widths[insertion_idx], nothing, (small_to_large, true, prev_res_idx, prev_res_idx))
                    transition_insertion_indices[b] = insertion_idx
                    b+=1
                    transition_data[b] = (spline_points[:, insertion_idx], q[:, insertion_idx], r[:, insertion_idx], s[:, insertion_idx], rectangle_widths[insertion_idx], nothing,  (small_to_large, false, res_idx, prev_res_idx))
                    transition_insertion_indices[b] = insertion_idx+1
                    b+=1                   
                end
                push!(fixed_indices, insertion_idx)

                if(curr_ss == BiochemicalAlgorithms.SecondaryStructure.SHEET) # skip array elements for arrow frames
                    arrow_insert_indices[a] = b
                    a+=1
                    b+=1
                end
            
            end
            prev_res_idx = res_idx
        end


        # arrows part 2 (add transition_point at the begin of the arrow head)
        a = 1
        for i=eachindex(arrow_starts)
            insertion_idx = arrow_starts[i]
            transition_data[arrow_insert_indices[a]] = (spline_points[:, insertion_idx], q[:, insertion_idx], r[:, insertion_idx], s[:, insertion_idx], T(1.0), sample_to_residue_indices[insertion_idx], nothing)
            transition_insertion_indices[arrow_insert_indices[a]] = insertion_idx+(n_to_c ? 0 : 1)
            a += 1
        end


    end



    # filter to decrease the amount of geometry to generate
    if(config.filter==Filter.ANGLE)
        # create a list of points that are fixed and cannot be removed by the filter
        push!(fixed_indices, 1, size(spline_points, 2)) # begin and end cannot be dropped

        if(config.color==Color.SECONDARY_STRUCTURE && config.backbone_type!=BackboneType.CARTOON)
            prev_res_idx = sample_to_residue_indices[1]
            prev_ss = spline.residue_info_dict[sample_to_residue_indices[1]][2]
            for (i, res_idx) in enumerate(sample_to_residue_indices)

                if(res_idx!=prev_res_idx && spline.residue_info_dict[res_idx][2]!=prev_ss)
                    push!(fixed_indices, i-1, i) # Color changes at secondary structure changes. 
                    prev_ss = spline.residue_info_dict[res_idx][2]
                    prev_res_idx = res_idx
                end
            end
        end
        if(config.color==Color.RESIDUE)
            prev_res_idx = sample_to_residue_indices[1]
            prev_ss = spline.residue_info_dict[sample_to_residue_indices[1]][2]

            prev_was_nothing = false
            for (i, resIdx) in enumerate(sample_to_residue_indices)

                if(resIdx===nothing)
                    prev_was_nothing = true
                    continue
                end


                if(resIdx!=prev_res_idx && !prev_was_nothing)
                    push!(fixed_indices, i-1, i) # Color changes at residue boundaries. The color should not be interpolated over a larger distance.
                    prev_res_idx = resIdx
                end

                prev_was_nothing = false
            end
        end

        local remaining_indices::Vector{Int}
        remaining_indices, remaining_count = filter_points_threshold(q, r, fixed_indices, with_color=(config.color==Color.RAINBOW))
    end

    log_info(types, "Type of spline points: ", typeof(spline_points))

    num_vertices = (remaining_count+num_transition_points)*config.resolution + 2 # end "caps"
    num_faces = (remaining_count+num_transition_points-1)*config.resolution*2 + 2*config.resolution # connections between neighbor circles + connections to the end points
    spline_mesh = PlainMesh(Array{T}(undef, 3, num_vertices), Array{T}(undef, 3, num_vertices), Array{Int}(undef, 3, num_faces), Vector{NTuple{3, Int}}(undef, num_vertices))
    # iterate and create vertices
    index_regular_frames = 1
    index_inserted_frames = 1
    i = 1 # currently highest index of index_regular_frames and the values of index_inserted_frames
    j = 1 # count of inserted frames (including the current one)
    while index_inserted_frames<=num_transition_points || index_regular_frames<=length(remaining_indices)
        if(config.color==Color.RAINBOW)
            fixed_color = rainbow(j/(remaining_count+num_transition_points))
        end

        if index_inserted_frames<=num_transition_points && transition_insertion_indices[index_inserted_frames] == i
            # insert this first
            
            # TODO resolution und filter koppeln
            @views generate_geometry_at_point!(spline_mesh, j,
                transition_data[index_inserted_frames][1],
                transition_data[index_inserted_frames][2],
                transition_data[index_inserted_frames][3],
                transition_data[index_inserted_frames][4], 
                spline.residue_info_dict,
                transition_data[index_inserted_frames][6],
                transition_data[index_inserted_frames][7],
                transition_data[index_inserted_frames][5],
                fixed_color,
                config)

            index_inserted_frames+=1
            j += 1
            continue
        end

        if (config.filter==Filter.NONE || (config.filter==Filter.ANGLE && remaining_indices[index_regular_frames]!=-1))
            # insert regular frame

            # TODO resolution und filter koppeln
            @views generate_geometry_at_point!(spline_mesh, j,
                spline_points[:, index_regular_frames], 
                q[:, index_regular_frames],
                r[:, index_regular_frames], 
                s[:, index_regular_frames], 
                spline.residue_info_dict,
                sample_to_residue_indices[index_regular_frames],
                nothing,
                config.backbone_type==BackboneType.CARTOON ? rectangle_widths[index_regular_frames] : T(1.0),
                fixed_color,
                config)


            index_regular_frames += 1
            i += 1
            j += 1
            continue
        end

        index_regular_frames += 1
        i+=1


        # # sanity check: frame should be orthogonal
        # @views if(!approx_zero(dot(q[:, current_index], r[:, current_index])) || !approx_zero(dot(q[:, current_index], s[:, current_index])) || !approx_zero(dot(s[:, current_index], r[:, current_index])))
        #     log_info(frame_rotation, current_index, " wrong angles ", dot(q[:, current_index], r[:, current_index]), " ", dot(q[:, current_index], s[:, current_index]), " ", dot(s[:, current_index], r[:, current_index]), " # ", q[:, current_index], " ", r[:, current_index], " ", s[:, current_index])
        # end
    end


    # add start and end vertices
    for i=axes(spline_points, 2)
        if(config.filter==Filter.NONE || remaining_indices[i]!=-1)
            
            spline_mesh.vertices[:, end-1] = spline_points[:,i]
            spline_mesh.normals[:, end-1] = -q[:, i]
            spline_mesh.colors[end-1] = spline_mesh.colors[1]
            break
        end
    end
    for i=reverse(axes(spline_points, 2))
        if(config.filter==Filter.NONE || remaining_indices[i]!=-1)
            spline_mesh.vertices[:, end] = spline_points[:,i]
            spline_mesh.normals[:, end] = q[:, i]
            spline_mesh.colors[end] = spline_mesh.colors[end-2]
            break
        end
    end


    add_faces_to_tube_mesh!(spline_mesh, config.resolution, remaining_count+num_transition_points)


    log_info(types, "Type of spline mesh: ", typeof(spline_mesh))

    # ----- debug export -----
    # cs = merge_multiple_meshes(circle_meshes)
    # export_mesh_representation_to_ply("circles.ply", Representation(cs))

    # fs = reduce(merge, framesA)
    # export_mesh_to_ply("framesA.ply", fs)

    # fs = reduce(merge, framesB)
    # export_mesh_to_ply("framesB.ply", fs)
    # ------------------------
    
    return spline_mesh
end

"""
Generates a PlainMesh for a system. 
"""
function prepare_backbone_model(
    ac::System{T}, config::BackboneConfig{T}) where {T<:Real}

    start_time = now()


    log_info(types, "Type: ", T)

    if(config.color==Color.ELEMENT)
        throw(ArgumentError("backbone-based models cannot be colored by elements of individual atoms"))
    end
    if(config.control_point_strategy==ControlPoints.C_ALPHA && config.frame==Frame.SECOND_SPLINE)
        throw(ArgumentError("for a second spline, ControlPoints.MID_POINTS is mandatory"))
    end

    

    if(config.color==Color.UNIFORM)
        uniform_color = (255, 0, 0)
    elseif(config.color==Color.CHAIN)
        chain_colors = n_colors(nchains(ac))
    end

    empty_mesh = PlainMesh(zeros(T, (3, 0)), zeros(T, (3, 0)), zeros(Int, (3, 0)), Vector{NTuple{3, Int}}())
    chain_meshes = fill(empty_mesh, (nchains(ac)))
    for (chain_idx, chain) in enumerate(BiochemicalAlgorithms.chains(ac))
        try
            color = nothing
            if(config.color==Color.UNIFORM)
                color = uniform_color
            elseif(config.color==Color.CHAIN)
                color = chain_colors[chain_idx]
            end

            chain_mesh = prepare_backbone_model(chain, config, color)

            chain_meshes[chain_idx] = chain_mesh
        catch e
            if e isa ErrorException
                log_warning("Skipped chain $(chain.name), because an error occured: $(e.msg)")
            else
                rethrow(e)
            end
        end
    end
    if(length(chain_meshes)==0)
        throw(ErrorException("No chain meshes were generated. "))
    end
    temp = merge_multiple_meshes(chain_meshes)
    log_info(types, "Type of result: ", typeof(temp))

    log_info(time_info, "Generated backbone mesh in $((now()-start_time).value/1000) seconds. ($(size(temp.vertices, 2)) vertices)")

    return temp

end