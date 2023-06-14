function shift!(mesh::Mesh, vector)
    for i=1:length(mesh.vertices)
        mesh.vertices[i] += vector
    end
end