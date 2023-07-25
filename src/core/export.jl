function export_mesh_representation_to_ply(path::AbstractString, representation::Representation)

    @assert isempty(representation.primitives) # Only representations containing a mesh are exported here. 
    @assert length(representation.vertices)%3 == 0

    colors = first(representation.colors).second
    @assert length(representation.vertices) == length(colors)*3

    stream = open(path, "w")

    #header
    println(stream, "ply")
    println(stream, "format ascii 1.0")
    println(stream, "element vertex $(length(representation.vertices) ÷ 3)")
    for dim in ["x", "y", "z"]
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
        r = parse(Int, colors[i][2:3], base=16)
        g = parse(Int, colors[i][4:5], base=16)
        b = parse(Int, colors[i][6:7], base=16)
        println(stream, "$x $y $z $r $g $b")
    end

    #Face List
    for i=1:(length(representation.connections) ÷ 3)
        println(stream, "3 $(representation.connections[(i-1)*3+1]) $(representation.connections[(i-1)*3+2]) $(representation.connections[(i-1)*3+3])")
    end

    close(stream)
end
