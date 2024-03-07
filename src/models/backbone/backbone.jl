export prepare_backbone_model

"""
Inserts elem into array. Assuming that array is sorted ascending, it will be still sorted after the insertion. 
"""
function insert_sorted!(array, elem)
    index = searchsortedfirst(array, elem)
    insert!(array, index, elem)
end

function check_config(user_config::Union{Nothing, PartialBackboneConfig}, T)
    default_config = BackboneConfig{T}(
        T(0.2), 
        T(1.5), 
        12, 
        BackboneType.BACKBONE, 
        Color.RAINBOW, 
        Spline.CUBIC_B, 
        ControlPoints.MID_POINTS, 
        Frame.RMF, 
        Filter.ANGLE)

    if user_config===nothing
        config = default_config
    else
        config = complete_config(user_config, default_config)
    end
    if(config.stick_radius<=0)
        throw(ArgumentError("stick_radius has to be >0"))
    end
    if(config.resolution_along<0.7)
        throw(ArgumentError("resolution along the spline has to be >=0.7")) # TODO why does it crash with 0.4, for example?
    end
    if(config.resolution_cross<3)
        throw(ArgumentError("at least three vertices per cross-section are necessary (resolution_cross>=3)"))
    end
    if(config.color==Color.ELEMENT)
        throw(ArgumentError("backbone-based models cannot be colored by elements of individual atoms"))
    end
    if(config.control_point_strategy==ControlPoints.C_ALPHA && config.frame==Frame.SECOND_SPLINE)
        throw(ArgumentError("for a second spline, ControlPoints.MID_POINTS is mandatory"))
    end

    return config
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
            if(frames_in_residue_count >= 2
                && (prev_res_idx ∈ arrow_fragment_indices || (is_last_frame && res_idx ∈ arrow_fragment_indices)))
                num_arrow = max(2, Int(round(frames_in_residue_count*2/3)))
                num_uniform = frames_in_residue_count - num_arrow
                
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
- `normal::AbstractVector{T}`
- `binormal::AbstractVector{T}`
- `linked_residue_idx::Int`: used to look up the residue in residue_info_dict (for secondary structure and color information)
- `rectangle_width::T`: Normally 1.0, except for frames that are part of an arrow head
- `fixed_color::Union{NTuple{3, Int}, Nothing}`: contains the fixed color if there is one (e. g. fixed color for the whole chain), otherwise nothing
- `residue_info_dict::Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}`: Maps indices of residues to tuples of (3 letter amino acid name, secondary structure)
- `config::BackboneConfig{T}`
"""
function generate_geometry_at_point!(
    result_mesh::PlainMesh{T},
    result_mesh_index::Int, 

    point::AbstractVector{T}, 
    normal::AbstractVector{T}, 
    binormal::AbstractVector{T}, 
    linked_residue_idx::Int, 
    rectangle_width::T, 

    fixed_color::Union{NTuple{3, Int}, Nothing}, 
    residue_info_dict::Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}},
    config::BackboneConfig{T}) where T
    # generate cross-section vertices
    start_index = (result_mesh_index-1)*config.resolution_cross+1
    end_index = result_mesh_index*config.resolution_cross
    if(config.backbone_type==BackboneType.BACKBONE)
        create_circle_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution_cross, config.stick_radius)
    elseif(config.backbone_type==BackboneType.RIBBON)
        create_ellipse_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution_cross, T(3)*config.stick_radius, config.stick_radius)
    elseif(config.backbone_type==BackboneType.CARTOON)
        structure = residue_info_dict[linked_residue_idx][2]
        if(structure==BiochemicalAlgorithms.SecondaryStructure.NONE)
            create_circle_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution_cross, config.stick_radius)
        elseif(structure==BiochemicalAlgorithms.SecondaryStructure.HELIX)
            create_ellipse_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution_cross, T(3)*config.stick_radius, T(1.5)*config.stick_radius)
        elseif(structure==BiochemicalAlgorithms.SecondaryStructure.SHEET)
            create_rectangle_in_local_frame!(@view(result_mesh.vertices[:, start_index:end_index]), @view(result_mesh.normals[:, start_index:end_index]), point, normal, binormal, config.resolution_cross, T(3)*config.stick_radius * rectangle_width, T(0.5)*config.stick_radius)
        end
    end

    # color
    color = nothing
    if(fixed_color!==nothing)
        color = fixed_color
    elseif(config.color==Color.SECONDARY_STRUCTURE)
        structure = residue_info_dict[linked_residue_idx][2]
        color = SS_COLORS[structure]
    elseif(config.color==Color.RESIDUE)
        aa = residue_info_dict[linked_residue_idx][1]
        color = AA_COLORS[aa]
    end
    for i=start_index:end_index
        result_mesh.colors[i] = color
    end
end

"""
Generates a single colored point and stores the vertices in result_mesh. 

# Arguments
- `result_mesh::PlainMesh{T}`: A preallocated struct that will contain the result. 
- `result_mesh_index::Int`: The created frame is the result_mesh_index-th frame along the spline. This parameter determines where in result_mesh's arrays, data is inserted. 
- `point::AbstractVector{T}`: The position of the frame
- `tangent::AbstractVector{T}`
- `linked_residue_idx::Int`: used to look up the residue in residue_info_dict (for secondary structure and color information)
- `fixed_color::Union{NTuple{3, Int}, Nothing}`: contains the fixed color if there is one (e. g. fixed color for the whole chain), otherwise nothing
- `residue_info_dict::Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}`: Maps indices of residues to tuples of (3 letter amino acid name, secondary structure)
- `config::BackboneConfig{T}`
"""
function generate_geometry_at_point!(
    result_mesh::PlainMesh{T},
    result_mesh_index::Int, 

    point::AbstractVector{T}, 
    tangent::AbstractVector{T}, 
    linked_residue_idx::Int,

    fixed_color::Union{NTuple{3, Int}, Nothing}, 
    residue_info_dict::Dict{Int, Tuple{String, BiochemicalAlgorithms.SecondaryStructure.T}},
    config::BackboneConfig{T}) where T

    start_index = (result_mesh_index-1)*config.resolution_cross+1
    end_index = result_mesh_index*config.resolution_cross

    # geometry
    result_mesh.vertices[:, start_index:end_index] = stack(repeat([point], config.resolution_cross)) # only 1 instead of resolution vertices would suffice, but then connect_circles_to_tube has to be modified (TODO)
    result_mesh.normals[:, start_index:end_index] = stack(repeat([tangent], config.resolution_cross))

    # color
    color = nothing
    if(fixed_color!==nothing)
        color = fixed_color
    elseif(config.color==Color.SECONDARY_STRUCTURE)
        structure = residue_info_dict[linked_residue_idx][2]
        color = SS_COLORS[structure]
    elseif(config.color==Color.RESIDUE)
        aa = residue_info_dict[linked_residue_idx][1]
        color = AA_COLORS[aa]
    end
    for i=start_index:end_index
        result_mesh.colors[i] = color
    end
end

"""
Generates a PlainMesh for chain. 
When the whole mesh should have a uniform color, it can be passed as fixed_color. 
"""
function prepare_backbone_model(chain::Chain{T}, partial_config::Union{PartialBackboneConfig, Nothing}=nothing; fixed_color::Union{Nothing, NTuple{3, Int}} = nothing) where {T<:Real}

    config = check_config(partial_config, T)

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
    elseif(config.spline==Spline.LINEAR)
        spline = Linear(chain, config.control_point_strategy)
    end

    # sample along spline
    spline_points, sample_to_residue_indices::Vector{Union{Int, Nothing}} = calculate_points(spline, config.resolution_along) #alloc
    velocities = calculate_velocities(spline, config.resolution_along)
    
    # construct local frames
    local q::Matrix{T}
    local r::Matrix{T}
    local s::Matrix{T}
    if(config.frame==Frame.RMF)
        q, r, s = rmf(spline_points, velocities)
    elseif(config.frame==Frame.SECOND_SPLINE)
        second_spline_points = calculate_minor_points(spline, config.resolution_along)
        q, r, s = frames_from_two_splines(spline_points, velocities, second_spline_points)
    end


    # when secondary structure is displayed, additional points are necessary
    ss_count = get_ss_count(sample_to_residue_indices, spline.residue_info_dict)
    num_transition_points = 0
    if(config.backbone_type==BackboneType.CARTOON)
        num_transition_points += 3*(ss_count[BiochemicalAlgorithms.SecondaryStructure.NONE]+ss_count[BiochemicalAlgorithms.SecondaryStructure.HELIX]+ss_count[BiochemicalAlgorithms.SecondaryStructure.SHEET]-1) # changes in secondary structure
        num_transition_points += ss_count[BiochemicalAlgorithms.SecondaryStructure.SHEET] # start of arrow heads
    end

    # The additional points will not be inserted into spline_points because of efficient memory allocation. 
    # Instead we store them (with related additional data) in a separate vector
    # The tuples contain (position, q, r, s, rectangle_width, residue_index, (is_small_to_large_transition, first_transition_part, residue_for_structure, residue_for_color))
    transition_data = Vector{Union{Tuple{Vector{T}, Vector{T}, Int}, Tuple{Vector{T}, Vector{T}, Vector{T}, Int, T}}}(undef, num_transition_points)
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
                if(n_to_c) # transition happens on the left side of the gap, because arrowheads cannot be strechted to the right

                    transition_data[b] = (spline_points[:, i-1], q[:, i-1], prev_res_idx) # the previous color in a single center-point
                    b+=1
                    transition_data[b] = (spline_points[:, i-1], q[:, i-1], res_idx) # the next color in a single center-point
                    b+=1
                    transition_data[b] = (spline_points[:, i-1], r[:, i-1], s[:, i-1], res_idx, T(1.0)) # a cross-section of the following type
                    b+=1
                else # transition happens on the right side of the gap
                    transition_data[b] = (spline_points[:, i], r[:, i], s[:, i], prev_res_idx, T(1.0)) # a cross-section of the previous type
                    b+=1
                    transition_data[b] = (spline_points[:, i], q[:, i], prev_res_idx) # the previous color in a single center-point
                    b+=1
                    transition_data[b] = (spline_points[:, i], q[:, i], res_idx) # the next color in a single center-point
                    b+=1
                end
                transition_insertion_indices[b-3:b-1] .= i

                push!(fixed_indices, i-1)
                push!(fixed_indices, i)

                

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
            transition_data[arrow_insert_indices[a]] = (spline_points[:, insertion_idx], r[:, insertion_idx], s[:, insertion_idx], sample_to_residue_indices[insertion_idx], T(1.0))
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

        remaining_indices, remaining_count = filter_points_threshold(q, r, fixed_indices, with_color=(config.color==Color.RAINBOW))
    else
        remaining_indices = 1:size(q, 2)
        remaining_count = size(q, 2)
    end

    log_info(types, "Type of spline points: ", typeof(spline_points))


    # allocate memory
    num_vertices = (remaining_count+num_transition_points)*config.resolution_cross + 2 # end "caps"
    spline_mesh = PlainMesh(Array{T}(undef, 3, num_vertices), Array{T}(undef, 3, num_vertices), Array{Int}(undef, 3, 0), Vector{NTuple{3, Int}}(undef, num_vertices))

    # weave spline_points and transition_points together
    picker = Array{Int}(undef, 2, remaining_count+num_transition_points)

    index_spline_points = 1
    index_transition_points = 1
    i = 1 # currently highest index of index_spline_points and the values of index_transition_points
    j = 1 # running variable
    while index_transition_points<=num_transition_points || index_spline_points<=length(remaining_indices)
        if index_transition_points<=num_transition_points && transition_insertion_indices[index_transition_points] == i
            picker[1, j] = 1
            picker[2, j] = index_transition_points

            index_transition_points+=1
            j += 1
            continue
        end


        if (config.filter==Filter.NONE || (config.filter==Filter.ANGLE && remaining_indices[index_spline_points]!=-1))
            picker[1, j] = 0
            picker[2, j] = index_spline_points

            index_spline_points += 1
            i += 1
            j += 1
            continue
        end

        index_spline_points += 1
        i+=1
    end

    # iterate and create vertices
    #
    Threads.@threads for j=1:remaining_count+num_transition_points
        local color_in_thread
        if(config.color==Color.RAINBOW)
            color_in_thread = rainbow(j/(remaining_count+num_transition_points))
        else
            color_in_thread = fixed_color
        end

        if picker[1, j]==0 # ordinary point
            # TODO resolution und filter koppeln
            @views generate_geometry_at_point!(spline_mesh, j,
                spline_points[:, picker[2, j]], 
                r[:, picker[2, j]], 
                s[:, picker[2, j]], 
                sample_to_residue_indices[picker[2, j]],
                config.backbone_type==BackboneType.CARTOON ? rectangle_widths[picker[2, j]] : T(1.0),
                color_in_thread,
                spline.residue_info_dict,
                config)
        else # transition point
            # TODO resolution und filter koppeln
            @views generate_geometry_at_point!(spline_mesh, j,
                transition_data[picker[2, j]]...,
                color_in_thread,
                spline.residue_info_dict,
                config)
        end

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


    add_faces_to_tube_mesh!(spline_mesh, config.resolution_cross, remaining_count+num_transition_points)


    log_info(types, "Type of spline mesh: ", typeof(spline_mesh))

    # ----- debug export -----
    # cs = merge_meshes(circle_meshes)
    # export_mesh_representation_to_ply("circles.ply", Representation(cs))

    # fs = reduce(merge, framesA)
    # export_mesh_to_ply("framesA.ply", fs)

    # fs = reduce(merge, framesB)
    # export_mesh_to_ply("framesB.ply", fs)
    # ------------------------
    
    return Representation(spline_mesh)
end

"""
Generates a PlainMesh for a system. 
"""
function prepare_backbone_model(ac::System{T}, partial_config::Union{PartialBackboneConfig, Nothing}=nothing) where {T<:Real}

    start_time = now()


    log_info(types, "Type: ", T)

    if(partial_config.color==Color.CHAIN)
        chain_colors = n_colors(nchains(ac))
    end

    chain_reps::Vector{Union{Missing, Representation{T}}} = fill(missing, (nchains(ac)))
    for (chain_idx, chain) in enumerate(BiochemicalAlgorithms.chains(ac))
        try
            color = nothing
            if(partial_config.color==Color.CHAIN)
                color = chain_colors[chain_idx]
            end

            chain_rep = prepare_backbone_model(chain, partial_config, fixed_color=color)

            chain_reps[chain_idx] = chain_rep
        catch e
            if e isa ErrorException
                log_warning("Skipped chain $(chain.name), because an error occured: $(e.msg)")
            else
                rethrow(e)
            end
        end
    end
    result_reps::Vector{Representation{T}} = filter(!ismissing, chain_reps)
    if(length(result_reps)==0)
        log_info(misc, "No chain meshes were generated. ")
    end
    result = merge_representations(result_reps)
    log_info(types, "Type of result: ", typeof(result))

    log_info(time_info, "Generated backbone mesh in $((now()-start_time).value/1000) seconds. ($(size(result.vertices, 2)) vertices)")

    return result

end


# convenience functions with default settings
function prepare_ribbon_model(
    ac::System{T}, partial_config::Union{PartialBackboneConfig, Nothing}=nothing) where {T<:Real}

    default = PartialBackboneConfig(0.2, 
        1.5, 12, 
        BackboneType.RIBBON, 
        Color.CHAIN, 
        Spline.CUBIC_B, 
        ControlPoints.MID_POINTS, 
        Frame.SECOND_SPLINE, 
        Filter.ANGLE)
    if partial_config!==nothing
        add_to_config!(partial_config, default)
    else
        partial_config = default
    end
    
    return prepare_backbone_model(ac, partial_config)
end

function prepare_cartoon_model(
    ac::System{T}, partial_config::Union{PartialBackboneConfig, Nothing}=nothing) where {T<:Real}

	default = PartialBackboneConfig(0.2, 
        1.5, 12, 
        BackboneType.CARTOON, 
        Color.SECONDARY_STRUCTURE, 
        Spline.CUBIC_B, 
        ControlPoints.MID_POINTS, 
        Frame.SECOND_SPLINE, 
        Filter.ANGLE)
    if partial_config!==nothing
        add_to_config!(partial_config, default)
    else
        partial_config = default
    end

    return prepare_backbone_model(ac, partial_config)
end