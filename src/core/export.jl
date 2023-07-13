function write_mesh_as_ply(path::AbstractString, mesh::Union{SimpleMesh, ColoredMesh})

    the_vertices = collect(mesh.vertices)
    faces = []
    for f in elements(topology(mesh))
        push!(faces, f.indices)
    end

    stream = open(path, "w")

    #header
    println(stream, "ply")
    println(stream, "format ascii 1.0")
    println(stream, "element vertex $(length(the_vertices))")
    for dim in ["x", "y", "z"]
        println(stream, "property float $dim") # TODO könnte auch double sein (4 oder 8 bytes)
    end
    
    if(typeof(mesh)<:ColoredMesh)
        for channel in ["red", "green", "blue"]
            println(stream, "property uchar $channel")
        end
    end

    println(stream, "element face $(length(faces))")
    println(stream, "property list uchar int vertex_index") # TODO correct datatypes
    println(stream, "end_header")

    #Vertex List
    if(typeof(mesh)<:ColoredMesh)
        for (v, c) in zip(the_vertices, mesh.colors)
            println(stream, "$(v.coords[1]) $(v.coords[2]) $(v.coords[3]) $(c[1]) $(c[2]) $(c[3])")
        end
    else
        for v in the_vertices
            println(stream, "$(v.coords[1]) $(v.coords[2]) $(v.coords[3])")
        end
    end

    #Face List
    for f in faces
        print(stream, "$(length(f))")
        for i in f
            print(stream, " $(i-1)")
        end
        println(stream)
    end

    close(stream)
end

function export_into_arrays(path::AbstractString, mesh::Union{SimpleMesh, ColoredMesh}, c_alpha_positions)
    

    num_point_coords = length(mesh.vertices) * 3
    points = Array{Float64}(undef, num_point_coords)
    a = 0
    for v in mesh.vertices
        points[a+1] = v.coords[1]
        points[a+2] = v.coords[2]
        points[a+3] = v.coords[3]
        
        a+=3
    end

    num_connects = nelements(topology(mesh)) * 3
    connections = Array{Int64}(undef, num_connects)  
    a = 0
    for f in elements(topology(mesh))
        if(length(f.indices)!=3)
            println("quad instead of tri", f, " ", a)
        else
        connections[a+1:a+3] = [convert.(Int64, (f.indices .-1))...]
        end
        a+=3
    end

    file = open(path, "w")
    print(file, "c_alpha_positions=[")
    for v in c_alpha_positions
        println(file, "($(v[1]), $(v[2]), $(v[3])), ")
    end
    println(file, "]")

    print(file, "vertices=")
    println(file, points)

    print(file, "connections=")
    println(file, connections)
    close(file)
end