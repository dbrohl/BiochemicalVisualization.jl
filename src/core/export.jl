export export_to_ply

"Calls an appropriate specific export functions for different geometry/mesh/representation types. "
function export_to_ply(path::AbstractString, representationOrMesh::Union{Meshes.SimpleMesh, ColoredMesh, Representation}; resolution=8)
    if typeof(representationOrMesh) <: Representation
        if ismesh(representationOrMesh)
            export_mesh_representation_to_ply(path, representationOrMesh)
        elseif isprimitivecollection(representationOrMesh)
            export_primitive_representation_to_ply(path, representationOrMesh, resolution)
        else
            throw(ArgumentError("Representation is neither mesh nor collection of primitives"))
        end
    else
        export_mesh_to_ply(path, representationOrMesh)
    end
end


function export_mesh_representation_to_ply(path::AbstractString, representation::Representation)
 
    @assert length(representation.vertices)%3 == 0

    colors = first(representation.colors).second
    @assert length(representation.vertices) == length(colors)*3

    stream = open(path, "w")

    #header
    println(stream, "ply")
    println(stream, "format ascii 1.0")
    println(stream, "element vertex $(length(representation.vertices) ÷ 3)")
    for dim in ["x", "y", "z", "nx", "ny", "nz"]
        println(stream, "property float $dim") # TODO könnte auch double sein (4 oder 8 bytes)
    end
    
    for channel in ["red", "green", "blue"]
        println(stream, "property uchar $channel")
    end

    println(stream, "element face $(length(representation.connections) ÷ 3)")
    println(stream, "property list uchar int vertex_index") # TODO correct datatypes
    println(stream, "end_header")

    #Vertex List
    for i=1:length(colors)
        x = representation.vertices[(i-1)*3 + 1]
        y = representation.vertices[(i-1)*3 + 2]
        z = representation.vertices[(i-1)*3 + 3]
        nx = representation.normals[(i-1)*3 + 1]
        ny = representation.normals[(i-1)*3 + 2]
        nz = representation.normals[(i-1)*3 + 3]
        r = parse(Int, colors[i][2:3], base=16)
        g = parse(Int, colors[i][4:5], base=16)
        b = parse(Int, colors[i][6:7], base=16)
        println(stream, x, " ", y, " ", z, " ", nx, " ", ny, " ", nz, " ", r, " ", g, " ", b)
    end

    #Face List
    for i=1:(length(representation.connections) ÷ 3)
        println(stream, "3 ", " ", representation.connections[(i-1)*3+1] , " ", representation.connections[(i-1)*3+2] , " ", representation.connections[(i-1)*3+3])
    end

    close(stream)
end

"Exports previously used mesh types like Meshes.SimpleMesh and ColoredMesh. "
function export_mesh_to_ply(path::AbstractString, mesh::X) where {X<:Union{Meshes.SimpleMesh, ColoredMesh}}

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

"Exports a representation that contains geometric primitives. "
function export_primitive_representation_to_ply(path::AbstractString, representation::Representation{T}, resolution) where {T}

    # make mesh
    empty_mesh = PlainMesh{T}(zeros(T, (3, 0)), zeros(T, (3, 0)), zeros(Int, (3, 0)), Vector{NTuple{3, Int}}())
    meshes = [empty_mesh]
    for (key, primitive_list) in representation.primitives
        local color_list
        if in(key, keys(representation.colors))
            color_list = representation.colors[key]
            while(length(color_list))<length(primitive_list)
                push!(color_list, "#FF0000")
            end
        else
            color_list = repeat(["#FF0000"], length(primitive_list))
        end
        for (primitive, color) in zip(primitive_list, color_list)
            try
                mesh = primitive_to_mesh(primitive, resolution)
                color!(mesh, hex_to_rgb(color))
                push!(meshes, mesh)
            catch ex
                @warn "Primitive of type $(typeof(primitive)) could not be turned into a mesh: $ex"
                rethrow(ex)
            end
        end
    end
    merged_mesh = merge_meshes(meshes)

    # export mesh
    rep = Representation(merged_mesh)
    export_mesh_representation_to_ply(path, rep)
end

"Converts a GeometryPrimitive to a PlainMesh. "
function primitive_to_mesh(prim::GeometryBasics.GeometryPrimitive{N, T}, resolution) where {N, T}
    if N!=3
        throw(ArgumentError("primitive was $N-D object instead of 3-D"))
    end
    if typeof(prim) <: GeometryBasics.Pyramid
        throw(ArgumentError("cannot handle pyramid primitives yet"))
    end
    if typeof(prim) <: GeometryBasics.HyperRectangle
        vertices = stack(GeometryBasics.coordinates(prim))
        face_iterator = GeometryBasics.faces(prim)
    else
        if typeof(prim) <: GeometryBasics.Sphere
            resolution = resolution ÷ 2
        end
        vertices = stack(GeometryBasics.coordinates(prim, resolution))
        face_iterator = GeometryBasics.faces(prim, resolution)
    end

    triangles = zeros(Int, (3, 0))
    if(!isempty(face_iterator))
        if typeof(first(face_iterator)) <: GeometryBasics.TriangleFace
            triangles = stack(face_iterator)
        elseif typeof(first(face_iterator)) <: GeometryBasics.QuadFace
            triangles = Matrix{Int}(undef, 3, 2*length(face_iterator))
            for (index, quad) in enumerate(face_iterator)
                triangles[:, 2*(index-1)+1] = quad[1:3]
                triangles[:, 2*index] = [quad[1] quad[3] quad[4]]
            end
        else
            throw(ErrorException("unknown type of faces $(typeof(face_iterator))"))
        end
    end



    if typeof(prim) <: GeometryBasics.Sphere
        normals = stack(GeometryBasics.normals(prim, resolution))
    else
        normals = zeros(T, (3, size(vertices, 2))) #TODO generate correct normals
        normals[1, :] .= 1
    end


    return PlainMesh{T}(vertices, normals, triangles, Vector{NTuple{3, Int}}(undef, size(vertices, 2)))

end