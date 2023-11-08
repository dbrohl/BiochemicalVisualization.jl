# This is very old and has to be modified a lot!
export prepare_backbone_model_gpu

function my_norm(a, b, c)
    return sqrt(a^2 + b^2 + c^2)
end

@kernel function normal_kernel_4(spline_points, normals, normal_lengths, valid_normals)
    I = @index(Global, Linear)

    if(I+1 <= size(spline_points, 2))
        for i=1:3
            normals[i, I] = spline_points[i, I+1] .- spline_points[i, I] #TODO coalesced
        end
        n_l = my_norm(normals[1, I], normals[2, I], normals[3, I])
        if(n_l< 10^-5)
            valid_normals[I] = 0
        else
            valid_normals[I] = 1
        end
        normal_lengths[I] = n_l
    else
        normal_lengths[I] = 0
        valid_normals[I] = 0
    end
end

@kernel function position_kernel_94(spline_points, normals, normal_lengths, new_indices, vertices, connections, stick_radius, resolution, ::Val{N}) where {N}
    I = @index(Global, Linear)
    U = eltype(spline_points)

    if(I<=length(normal_lengths) && normal_lengths[I]>= 10^-5)
   
        center = @private U (3, )
        normal = @private U (3, )
        points = @private U (3, N)
        
        for i=1:3
            center[i] = spline_points[i, I]
            normal[i] = normals[i, I]
        end
        for i=1:resolution
            α = U(i/resolution*2*π)
            points[1, i] = stick_radius*cos(α)
            points[2, i] = stick_radius*sin(α)
            points[3, i] = 0
        end

        if(normal[1]!=0 && normal[2]!=0)
            if(normal[1]!=0)
                rotationAxis = Vec(-normal[2]/normal[1], U(1), U(0)) * -sign(normal[1])
            elseif(normal[2]!=0)
                rotationAxis = Vec(U(1), -normal[1]/normal[2], U(0)) * sign(normal[2]) # rotation axis is perpendicular to direction and lies in the xy-plane
            end

            rotationAngle = Meshes.∠(Vec{3,U}(normal[1], normal[2], normal[3]), Vec(U(0), U(0), U(1)))
            points = (points' * AngleAxis(rotationAngle, rotationAxis...))'
        end
        
        baseIndex = (new_indices[I]-1)*2*resolution+1
        baseVertexNum = (new_indices[I]-1)*resolution + 1
        r1 = @private Int (N, )
        r2 = @private Int (N, )

        r1 = UnitRange(baseVertexNum, baseVertexNum+resolution-1)
        r2 = UnitRange(baseVertexNum+resolution, baseVertexNum+2*resolution-1)

        for i=1:resolution
            vertices[(new_indices[I]-1)*resolution + i, 1] = center[1] + points[1, i]
            vertices[(new_indices[I]-1)*resolution + i, 2] = center[2] + points[2, i]
            vertices[(new_indices[I]-1)*resolution + i, 3] = center[3] + points[3, i]

            connections[baseIndex+i-1, 1] = r1[i]
            connections[baseIndex+i-1, 2] = r1[mod1(i+1, resolution)]
            connections[baseIndex+i-1, 3] = r2[i]

            connections[baseIndex+resolution+i-1, 1] = r2[i]
            connections[baseIndex+resolution+i-1, 2] = r2[mod1(i+1, resolution)]
            connections[baseIndex+resolution+i-1, 3] = r1[mod1(i+1, resolution)]
        end
    end

end

function prepare_backbone_model_gpu(
    ac::AbstractAtomContainer{T}; 
    stick_radius=T(0.2), resolution::Int=30) where {T<:Real}

    start_time = now()
    U = Float64
    if(T <: AbstractFloat)
        U = T
    end

    log_info(types, "Types: ", T, " ", U)

    vertices_per_unit = resolution / 2*π*stick_radius
    circle_prototype = discretize(Sphere(Point(U(0),U(0)), stick_radius), RegularDiscretization(resolution))
    circle_mesh = PlainNonStdMesh(lift_into_3d(circle_prototype))

    cap_mesh = PlainMesh(Hemisphere(stick_radius, resolution, U))

    log_info(types, "Types of proto meshes: ", typeof(circle_mesh), " ", typeof(cap_mesh))
    chain_meshes::Vector{PlainMesh{U}} = []
    chain_colors = map(c->map(channel->Int(channel*255), (c.r, c.g, c.b)), collect(distinguishable_colors(nchains(ac)+1))[2:end])


    for (chain_num, chain) in enumerate([collect(eachchain(ac))[1]])
        println("chain $chain_num")

        c_alphas = filter(x -> x.element==Elements.C && x.name=="CA", atoms(chain))
        @assert length(c_alphas)>=2 # TODO was sonst?

        c_alpha_spline = CatmullRom(hcat(map(x->x.r, c_alphas)...))
        
        spline_points::AbstractMatrix{T} =  c_alpha_spline(vertices_per_unit)
        # spline_points::AbstractMatrix{T} = [0 1 2 3 4
        #                                     0 0 0 0 0
        #                                     0 0 0 0 0]
        log_info(types, "Type of spline points: ", typeof(spline_points))
        

        spline_points_d = allocate(backend, U, size(spline_points)...)
        normals_d = allocate(backend, U, 3, size(spline_points, 2)-1)
        normal_lengths_d = allocate(backend, U, size(spline_points, 2)-1)
        valid_normals_d = allocate(backend, Int, size(spline_points, 2)-1)
        copyto!(spline_points_d, spline_points)

        println(size(spline_points), " ", size(spline_points_d), " ", size(valid_normals_d))

        k = normal_kernel_4(backend, 1024) # TODO try different values
        k(spline_points_d, normals_d, normal_lengths_d, valid_normals_d, ndrange=size(normals_d, 2))
        KernelAbstractions.synchronize(backend)

        
        new_indices_d = accumulate(+, valid_normals_d)
        last_val = CUDA.@allowscalar new_indices_d[end] #TODO platform-independant
        println("finished 1 ", last_val)
        vertices_d = allocate(backend, U, resolution*last_val, 3)
        connections_d = allocate(backend, Int, last_val*(2*resolution), 3)
        println("Size of connections: ", size(connections_d), " Size of spline_points: ", size(spline_points))
        k1 = position_kernel_94(backend, 1024)
        k1(spline_points_d, normals_d, normal_lengths_d, new_indices_d, vertices_d, connections_d, stick_radius, resolution, Val(resolution),  ndrange=size(spline_points, 2))
        KernelAbstractions.synchronize(backend)

        export_vertices = convert(Array, vertices_d)'
        export_connects = convert(Array, connections_d)'
        println(size(export_vertices))
        temp = [Point3(v...) for v in eachcol(export_vertices)]
        println("temp", size(temp))
        connects::Vector{Connectivity} = [connect(Tuple(x)) for x in eachcol(export_connects)]
        # connects::Vector{Connectivity} = [connect(Tuple((i-1)*(resolution+2)+1:(i-1)*(resolution+2)+resolution)) for i in 1:last_val]
        # connects::Vector{Connectivity} = [connect((1, 2))]
        # for i in 1:last_val
        #     push!(connects, connect((i*(resolution+2)-1, i*(resolution+2))))
        # end

        println("Connects:", connects[1:10])

        export_mesh_to_ply("gpu.ply", Meshes.SimpleMesh(temp, connects))

        # TODO caps, color, shift/flip

        # push!(chain_meshes, spline_mesh)
    end


    # temp = merge_multiple_meshes(chain_meshes)


    
    # result = Representation(temp)
    # log_info(types, "Type of result: ", typeof(result))

    # log_info(time_info, "Generated backbone mesh on gpu in $((now()-start_time).value/1000) seconds. ($(length(result.vertices) ÷ 3) vertices)")

    # return result
    return 0

end