function write_mesh_as_ply(path::AbstractString, mesh::SimpleMesh)::void
    the_vertices = collect(vertices(mesh))
    faces = []
    for f in elements(topology(mesh))
        push!(faces, f.indices)
    end
    println(faces)

    stream = open(path, "w")

    #header
    println(stream, "ply")
    println(stream, "format ascii 1.0")
    println(stream, "element vertex $(length(the_vertices))")
    for dim in ["x", "y", "z"]
        println(stream, "property float $dim") # TODO könnte auch double sein (4 oder 8 bytes)
    end
    println(stream, "element face $(length(faces))")
    println(stream, "property list uchar int vertex_index") # TODO correct datatypes
    println(stream, "end_header")

    #Vertex List
    for v in the_vertices
        println(stream, "$(v.coords[1]) $(v.coords[2]) $(v.coords[3])")
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