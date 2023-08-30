export gpu_rotation_test

@kernel function kernel_8(vertices, directions, outer_radius, num_points, ::Val{N}) where {N}
    I = @index(Global, Linear)
    U = Float32

    if(I<=size(directions, 2))

        #@print(I, " ", N, " ", num_points, "\n")
        normal = @private U (3, )
        points = @private U (3, N)

        for j=1:num_points
            for i=1:3
                points[i, j] = vertices[i, (I-1)*num_points+j]
            end
        end

        for i=1:3
            normal[i] = directions[i, I]
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

        for j=1:num_points
            for i=1:3
                vertices[i, (I-1)*num_points+j] = points[i, j] + normal[i]*outer_radius
            end
        end
    end

end

function gpu_rotation_test()
    stick_radius = 0.5
    resolution=15
    U = Float32

    circle_prototype = discretize(Sphere(Point(U(0),U(0)), stick_radius), RegularDiscretization(resolution))
    circle_mesh = PlainNonStdMesh(lift_into_3d(circle_prototype))
    all_points = [circle_mesh.vertices [0;0;0]]
    all_connects = [(1:resolution)'; circshift(1:resolution, 1)'; repeat([resolution+1], resolution)']
    all_colors = repeat([(100, 100, 100)], resolution+1)
    circle_std_mesh = PlainMesh(all_points, all_connects, all_colors)

    all_perm_meshes::Vector{PlainMesh{U}} = []
    permutations = [[1, 2, 3]]
    # [1, 3, 2],
    # [2, 1, 3],
    # [2, 3, 1],
    # [3, 1, 2],
    # [3, 2, 1]]

    for (height, perm) in enumerate(permutations)

        res = 16
        r = collect(range(0, 2*π, length=res+1))[1:end-1]
        directions = Matrix{U}(undef, 3, 0)
        colors = []

        for (i, z_angle) in enumerate(r)
            layer_radius = cos(z_angle)
            z = sin(z_angle)

            for (j, xy_angle) in enumerate(r)
                x = cos(xy_angle) * layer_radius
                y = sin(xy_angle) * layer_radius
                directions = [directions [x; y; z][perm]]
                push!(colors, convert.(Int, floor.(((x, y, z) .+ 1).*128)))
            end
        end


        vertices_d = allocate(backend, U, 3, size(circle_std_mesh.vertices, 2)*size(directions, 2))
        directions_d = allocate(backend, U, size(directions)...)

        copyto!(directions_d, directions)
        for i=1:size(directions, 2)
            copyto!(vertices_d, (i-1)*size(circle_std_mesh.vertices, 2)+1, circle_std_mesh.vertices, 1, length(circle_std_mesh.vertices))
        end

        #colors_d = allocate(backend, Int, 3, length(colors))
        println("parallel executions: ", size(directions, 2))
        k = kernel_8(backend, 256)
        k(vertices_d, directions_d, 5, size(circle_std_mesh.vertices, 2), Val(size(circle_std_mesh.vertices, 2)), ndrange=size(directions_d, 2))
        KernelAbstractions.synchronize(backend)

        sms::Vector{PlainMesh{U}} = []
        vertices_host = convert(Array, vertices_d)
        for i=1:length(colors)
            mesh = deepcopy(circle_std_mesh)
            mesh.vertices = vertices_host[:, (i-1)*resolution+1 : i*resolution]
            color!(mesh, colors[i])
            push!(sms, mesh)

            export_mesh_representation_to_ply("rotation_test_$i.ply", Representation(mesh))
        end

        temp = merge_multiple_meshes(sms)
        push!(all_perm_meshes, temp)

    end
    export_mesh_representation_to_ply("rotation_test.ply", Representation(merge_multiple_meshes(all_perm_meshes)))

end