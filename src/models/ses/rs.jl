# create reduced surface
    # create list of spheres(atom positions and radii+probe_radius (+ ggf. offset))
    # preprocessing
        # determine bounding box of atom centers (without vdw/probe radius)
        # determine max radius
        # create hash grid (Sanner et al use a tree instead of a grid. )

function create_reduced_surface(
    ac::System{T}, probe_radius::T) where {T<:Real}
    if natoms(ac)==0
        return ErrorException("No atoms in system. ")
    end

    # ----- void RSComputer::preProcessing() ------
    epsilon = 1e-4;

    first_atom = first(eachatom(ac))
    max_radius = 0
    bounding_box = [first_atom.r;; first_atom.r] # [x/y/z, min/max]

    for atom in eachatom(ac)
        if atom.radius == 0 # fix missing radius data once
            atom.radius = 1
        end
        max_radius = max(max_radius, atom.radius)

        bounding_box[:, 1] = min.(bounding_box[:, 1], atom.r)
        bounding_box[:, 2] = max.(bounding_box[:, 2], atom.r)
    end
    log_info(rs, "Bounding box of all $(natoms(ac)) atoms", bounding_box)

    max_dist = 2*(max_radius + probe_radius)
    hashgrid = HashGrid{T, AtomWrapper{T}}(bounding_box, [max_dist, max_dist, max_dist], atom -> atom.r, natoms(ac)) #TODO eventuell bounding_box anpassen und vergrößern (->reducedSurface.C, line 1861f.)
    # TODO ggf. mit CellListMap machen, falls das schneller ist

    # almost all atoms are inserted into the hashgrid. Exception: They are too close to an already inserted hypothetical neighbour
    for atom in eachatom(ac)
        too_close = false
        for (neighbour, _) in each_neighbour(hashgrid, atom.r)
            if (norm(atom.r .- neighbour.r) - 0.05*max(neighbour.radius, atom_radius)) < epsilon
                too_close = true
                break
            end
        end

        if !too_close
            push!(hashgrid, AtomWrapper(atom))
        end
    end
    log_info(rs, "Hashgrid contains $(length(hashgrid.data)) atoms")


    # cache actual neighbours
    for i = length(grid.data):-1:2
        offset = grid.data[i].atom.radius + 2*probe_radius
        pos = grid.data[i].r
        
        for (neighbour, dist) in each_neighbour(grid, grid.data[i].r)
            if neighbour.index_in_list > i # we only need to count every pair twice

                max_dist = (neighbour.atom.radius+offset)^2
                if !isgreater_tolerance(dist, max_dist)
                    push!(grid.data[i].neighbours, neighbour.index_in_list)
                    push!(neighbour.neighbours, i)
                end
            end
        end
    end

    # ----- reamaining part of void RSComputer::run() -----
    start = 1;
    while start!=0
        start = getStartPosition(grid, neighbourDict, probe_positions, probe_radius)

        if start==2
            extendComponent(new_vertices, rm_vertices, all_atoms, neighbourdict, probe_radius, new_faces, vertices, faces, number_of_faces, rs, probe_positions) #TODO resolve call
        elseif start==3
            getRSComponent(faces, number_of_faces, new_faces, all_atoms, neighbourdict, probe_radius, vertices, rs, rm_vertices, new_vertices, probe_positions) #TODO resolve call
        end
    end

    clean_rs(number_of_vertices, vertices, number_of_edges, edges, number_of_faces, faces) #TODO resolve function call

    #TODO if Constants::EPSILON is used later, it has to be the old value from line 676 in reducedSurface.C



end

function clean_rs(number_of_vertices, vertices, number_of_edges, edges, number_of_faces, faces)

    while ((number_of_vertices > 0) &&
        (vertices[number_of_vertices-1] === nothing))
        pop!(vertices)
        number_of_vertices-=1
    end

    for i=1:number_of_vertices
        if (vertices[i] === nothing)
            vertices[i] = vertices[number_of_vertices-1]
            vertices[i].index = i
            pop!(vertices)
            number_of_vertices -= 1
            while (vertices[number_of_vertices-1] === nothing)
                pop!(vertices)
                number_of_vertices -= 1
            end
        end
    end
    while ((number_of_edges > 0) &&
            (edges[number_of_edges-1] === nothing))
        pop!(edges)
        number_of_edges -= 1
    end
    for i=1:number_of_edges
        if (edges[i] === nothing)
            edges[i] = edges[number_of_edges_1]
            edges[i].index = i;
            edges_.pop_back();
            number_of_edges_--;
            while (edges_[number_of_edges_-1] == NULL)
                pop!(edges)
                number_of_edges -= 1
            end
        end
    end
    while ((number_of_faces > 0) &&
            (faces[number_of_faces-1] === nothing))
        pop!(faces)
        number_of_faces -= 1
    end
    for i=1:number_of_faces
    if (faces[i] === nothing)
        faces[i] = faces[number_of_faces-1]
        faces[i].index = i
        pop!(faces)
        number_of_faces -= 1
        while (faces[number_of_faces-1] === nothing)
            pop!(faces)
            number_of_faces -= 1
        end
    end
end



function getStartPosition(grid, neighbourDict, probe_positions, probe_radius, rs, new_vertices, vertices)
    if findFirstFace(grid, neighbourDict, probe_positions, probe_radius)!==nothing
        return 3
    end
    if findFirstEdge(grid, probe_radius, neighbourDict, rs, new_vertices, vertices)!==nothing
        return 2
    end
    if findFirstVertex(grid.data, rs, new_vertices, vertices)!==nothing
        return 1
    end
    return 0
end

function findFirstFace(grid, neighbourDict, probe_positions, probe_radius)
    for direction = 1:3
        for extreme = 1:2
            face = findFace(direction, extreme, grid, neighbourDict, probe_positions, probe_radius)
            if face!==nothing
                return face
            end
        end
    end
    return nothing
end
function findFirstEdge(grid, probe_radius, neighbourDict, rs, new_vertices, vertices)
    for direction = 1:3
        for extreme = 1:2
            edge = findEdge(direction, extreme, grid, probe_radius, neighbourDict, rs, new_vertices, vertices)
            if edge!==nothing
                return edge
            end
        end
    end
    return nothing
end
function findFirstVertex(all_atoms::AbstractVector{AtomWrapper}, rs, new_vertices, vertices)
    for (index, atom) in enumerate(all_atoms)
        if atom.status==UNKNOWN
            if length(atom.neighbours)==0
                vertex = RSVertex(index)
                insert(vertex, rs, new_vertices, vertices)
                all_atoms[index].status = AtomStatus.ON_SURFACE
                
                return vertex
            end
        end
    end
    return nothing
end

function findFace(direction, extreme, grid::HashGrid{AtomWrapper}, neighbourDict, probe_positions, probe_radius::T,rs, new_faces, new_vertices, vertices) where T
    idx1 = findFirstAtom(direction, extreme, grid.data)
    if idx1 === nothing
        return nothing
    end

    idx2 = findSecondAtom(idx1, direction, extreme, grid, probe_radius)
    if idx2 === nothing
        return nothing
    end

    possible_neighbours = neighboursOfTwoAtoms(grid.data[idx1], grid.data[idx2], neighbourDict)

    candidates = findThirdAtom(grid.data[idx1], grid.data[idx2], map(idx -> grid.data[idx], possible_neighbours), probe_positions, probe_radius)
    if length(candidates)==0
        return nothing
    end
    
    i = 1
    idx3 = -1
    found = false
    while (!found && i <= length(candidates))
        idx3, probe = candidates[i]
        found = ((grid.data[idx3] == AtomStatus.UNKNOWN)
                    && checkProbe(probe, grid.data[idx1], grid.data[idx2], grid.data[idx3], neighbourDict, grid.data))
        i += 1
    end
    if found
        vertex1 = RSVertex(idx1)
        vertex2 = RSVertex(idx2)
        vertex3 = RSVertex(idx3)
        
        e1 = RSEdge()
        e2 = RSEdge()
        e3 = RSEdge()
        
        face = RSFace{T}()
        
        connectRSElements!(vertex1,vertex2,vertex3,e1,e2,e3,face,probe, grid.data[idx1], grid.data[idx2], grid.data[idx3])
        
        insert(face, rs, new_faces)
        insert(vertex1, rs, new_vertices, vertices)
        insert(vertex2, rs, new_vertices, vertices)
        insert(vertex3, rs, new_vertices, vertices)
        grid.data[idx1].status = AtomStatus.ON_SURFACE
        grid.data[idx2].status = AtomStatus.ON_SURFACE
        grid.data[idx3].status = AtomStatus.ON_SURFACE

        
        return face;
    else
        grid.data[idx1] = AtomStatus.INSIDE
        grid.data[idx2] = AtomStatus.INSIDE
        return NULL;
    end

end

function findEdge(direction, extreme, grid, probe_radius, neighbourdict, rs, new_vertices, vertices)
    idx1 = findFirstAtom(direction,extreme, grid.data);
    if (idx1 == -1)
        return nothing
    end

    idx2 = findSecondAtom(idx1,direction,extreme, grid, probe_radius)
    if (idx2 == -1)
        return nothing
    end

    vertex1 = RSVertex(idx1)
    vertex2 = RSVertex(idx2)
    neighboursOfTwoAtoms(grid.data[idx1], grid.data[idx2], neighbourdict)

    edge = createFreeEdge(vertex1,vertex2, grid.data, probe_radius, neighbourdict)
    if (edge !== nothing)
        insert(edge, rs)
        insert(vertex1, new_vertices, vertices)
        insert(vertex2, new_vertices, vertices)
        grid.data[idx1].status = AtomStatus.ON_SURFACE
        grid.data[idx2].status = AtomStatus.ON_SURFACE

        return edge
    else
        delete!(grid.data[idx1].neighbours, idx2)
        delete!(grid.data[idx2].neighbours, idx1)

        return nothing
    end
end

function findFirstAtom(direction, extreme, all_atoms::AbstractVector{AtomWrapper})
    extreme_atom_index = nothing
    # find the first atom of unknown status
    index = nothing
    for (i, atom) in enumerate(all_atoms)
        if atom.status == UNKNOWN
            index = i
            break
        end
    end

    
    if index!==nothing
        extreme_atom_index = index
        next_atom = all_atoms[index]
        extreme_value = ((extreme == 1) ? next_atom.atom.r[direction]-next_atom.atom.radius 
                                        : next_atom.atom.r[direction]+next_atom.atom.radius)

        #find the atom of unknown status lying on the extrem position
        for i=index+1:length(all_atoms)
            if all_atoms[i].status == UNKNOWN
                next_atom = all_atoms[i]
                temp = ((extreme == 1) ? next_atom.atom.r[direction]-next_atom.atom.radius
                                       : next_atom.atom.r[direction]+next_atom.atom.radius)
                if (((extreme == 1) && isless_tolerance(temp, extreme_value)) ||
                        ((extreme != 1) && isgreater_tolerance(temp, extreme_value)))
                    extreme_value = temp
                    extreme_atom_index = i
                end
            end
        end
    end
    return extreme_atom_index
end

function findSecondAtom(atom_index_1, direction, extreme, grid, probe_radius)
    second_atom_index = nothing
    # find the first neighbour atom of unknown status
    selected_atom = nothing
    iter = each_neighbour(grid, grid.data[atom_index_1].atom.r)
    for (neighbour, distance) in iter
        if neighbour.status == UNKNOWN
            selected_atom = neighbour
            break
        end
    end



    if selected_atom!==nothing
        second_atom_index = selected_atom.index_in_list
        first_atom = grid.data[atom_index_1]
        first_atom.atom.radius += probe_radius
        extreme_value = ((extreme == 1) ? first_atom.atom.r[direction]+first_atom.atom.radius
                                        : first_atom.atom.r[direction]-first_atom.atom.radius)
        next_atom = nothing
        # find the neighbour atom of unknown status lying on the extreme position
        for (neighbour, distance) in iter
            if neighbour.status == UNKNOWN
                next_atom = neighbour
                next_atom.radius += probe_radius; # ensure that this doesnt destroy the underlying data
                intersection_circle = getIntersectionCircle(first_atom.atom ,next_atom.atom)
                if intersection_circle!==nothing
                    next_extreme = getCircleExtremum(intersection_circle,direction,extreme)
                    if (((extreme == 1) && isless_tolerance(next_extreme,extreme_value)) ||
                            ((extreme != 1) && isgreater_tolerance(next_extreme,extreme_value)))
                        extreme_value = next_extreme
                        second_atom_index = neighbour.index_in_list
                    end
                end
            end
        end
    end
    return second_atom_index
end

""" Returns the intersection circle or nothing. """
function getIntersectionCircle(a::Atom{T}, b::Atom{T}) where T
    norm = a.r .- b.r
    square_dist = dot(norm, norm)
    if (iszero_tolerance(square_dist))
        return nothing
    end
    dist = sqrt(square_dist)
    if isless_tolerance(a.radius + b.radius, dist)
        return nothing
    end
    if isgreaterorequal_tolerance(abs(a.radius - b.radius), dist)
        return nothing
    end

    radius1_square = a.radius^2
    radius2_square = b.radius^2
    u = radius1_square - radius2_square + square_dist
    length = u / (2 * square_dist)
    square_radius = radius1_square - u * length / 2
    if square_radius < 0
        return nothing
    end

    return Circle{T}(a.r + (norm * length), norm / dist, sqrt(square_radius))
end

"Returns a circle or nothing"
function getIntersectionCircle(sphere::ProbeSphere, plane::Plane)
    distance = getDistance(sphere.r, plane)
    if isgreater_tolerance(distance, sphere.radius)
        return nothing
    end

    v3 = normalize(plane.n)
    if isequal_tolerance(distance, sphere.radius)
        return Circle(sphere.r .+ sphere.radius.*v3, plane.n, 0)
    else
        return Circle(sphere.r .+ distance .* v3, plane.n,  sqrt(sphere.radius*2 - distance^2))
    end
end


"""Returns two intersection points or (nothing, nothing). """
function getIntersectionPoints(s1::Atom{T}, s2::Atom{T}, s3::Atom{T}, test::Bool=true)

    r1_square = s1.radius^2
    r2_square = s2.radius^2
    r3_square = s3.radius^2
    p1_square_length = dot(s1.r, s1.r)
    p2_square_length = dot(s2.r, s2.r)
    p3_square_length = dot(s3.r, s3.r)
    u = (r2_square-r1_square-p2_square_length+p1_square_length)/2
    v = (r3_square-r1_square-p3_square_length+p1_square_length)/2
    local plane1
    local plane2
    try
        plane1 = Plane{T}(s2.r .- s1.r, u)
        plane2 = Plane{T}(s3.r .- s1.r, v)
    catch e
        return nothing, nothing
    end

    line = getIntersectionLine(plane1,plane2)
    if line!==nothing
        diff .= s1.r .- line.r
        num_solutions, x1, x2 = SolveQuadraticEquation(dot(line.d, line.d), dot(-diff, line.d)*2, dot(diff, diff)-r1_square)
        if (num_solutions > 0)
            @. p1 = line.r + x1*line.d
            @. p2 = line.r + x2*line.d
            if test
                testVec .= s1.r .- p1
                if isnotequal_tolerance(test^2, r1_square)
                    return nothing, nothing
                end
                testVec .= s1.r .- p2
                if isnotequal_tolerance(test^2, r1_square)
                    return nothing, nothing
                end
                testVec .= s2.r .- p1
                if isnotequal_tolerance(test^2, r2_square)
                    return nothing, nothing
                end
                testVec .= s2.r .- p2
                if isnotequal_tolerance(test^2, r2_square)
                    return nothing, nothing
                end
                testVec .= s3.r .- p1
                if isnotequal_tolerance(test^2, r3_square)
                    return nothing, nothing
                end
                testVec .= s3.r .- p2
                if isnotequal_tolerance(test^2, r3_square)
                    return nothing, nothing
                end
            end
            return p1, p2
        end
    end
    return nothing, nothing
end

"""Returns two intersection points or (nothing, nothing). """
function getIntersectionPoints(probe::ProbeSphere, line::Line)
    number_of_solutions, x1, x2 = SolveQuadraticEquation(dot(line.d, line.d), dot((line.r .- probe.r), line.d .* 2), dot((line.r .- probe.r), (line.r .- probe.r)) - probe.radius^2)

    if (number_of_solutions == 0)
        return nothing, nothing
    end

    return line.r + x1*line.d, line.r + x2*line.d
end

"""
Computes the intersection line between plane1 and plane2. 
Returns either a Line or nothing. 
"""
function getIntersectionLine(plane1::Plane{T}, plane2::Plane{T})
    u = plane1.r*plane1.n
    v = plane2.r*plane2.n
    det = plane1.n[1]*plane2.n[2] - plane1.n[2]*plane2.n[1]
    line = Line()
    if iszero_tolerance(det)
        det = plane1.n[1]*plane2.n[3] - plane1.n[3]*plane2.n[1]
        if iszero_tolerance(det)
            det = plane1.n[2]*plane2.n[3] - plane1.n[3]*plane2.n[2]
            if iszero_tolerance(det)
                return nothing
            else
                a = plane2.n[3]/det
                b = -plane1.n[3]/det
                c = -plane2.n[2]/det
                d = plane1.n[2]/det
                line.r[1] = 0
                line.r[2] = a*u + b*v
                line.r[3] = c*u + d*v
                line.d[1] = -1
                line.d[2] = a*plane1.n[1] + b*plane2.n[1]
                line.d[3] = c*plane1.n[1] + d*plane2.n[1]
            end
        else
            a = plane2.n[3]/det
            b = -plane1.n[3]/det
            c = -plane2.n[1]/det
            d = plane1.n[1]/det
            line.r[1] = a*u + b*v
            line.r[2] = 0
            line.r[3] = c*u + d*v
            line.d[1] = a*plane1.n[2] + b*plane2.n[2]
            line.d[2] = -1
            line.d[3] = c*plane1.n[2] + d*plane2.n[2]
        end
    else
        a = plane2.n[2]/det
        b = -plane1.n[2]/det
        c = -plane2.n[1]/det
        d = plane1.n[1]/det
        line.r[1] = a*u + b*v
        line.r[2] = c*u + d*v
        line.r[2] = 0
        line.d[1] = a*plane1.n[2] + b*plane2.n[2]
        line.d[2] = c*plane1.n[2] + d*plane2.n[2]
        line.d[2] = -1
    end
    return line
end


function getCircleExtremum(circle::Circle{T}, direction, extreme)
    min = 0
    max = 0

    norm2 = circle.n .* circle.n

    if direction==0
        if iszero_tolerance(circle.n[2]) && iszero_tolerance(circle.n[3])
            min = max = circle.r[1]
        else
            x_norm = norm2[2] + norm2[3]
            x_norm /= norm2[1]+x_norm
            x_norm = circle.radius * sqrt(x_norm)
            min = (circle.r[1]) - x_norm
            max = (circle.r[1]) + x_norm
        end
    elseif direction==1
        if iszero_tolerance(circle.n[1]) && iszero_tolerance(circle.n[3])
            min = max = circle.r[2]
        else
            y_norm = norm2[1] + norm2[3]
            y_norm /= norm2[2] + y_norm
            y_norm = circle.radius * sqrt(y_norm)
            min = (circle.r[2])-y_norm
            max = (circle.r[2])+y_norm
        end
    elseif direction==2
        if iszero_tolerance(circle.n[1]) && iszero_tolerance(circle.n[2])
            min = max = circle.r[3]
        else
            z_norm = norm2[1] + norm2[2]
            z_norm /= norm2[3] + z_norm
            z_norm = circle.radius * sqrt(z_norm)
            min = circle.r[3] - z_norm
            max = circle.r[3] + z_norm
        end
    end
    if extreme == 0
        return min
    else
        return max
    end
end

function getDistance(point::AbstractVector{T}, plane::Plane{T}) where T
    len = norm(plane.n)
    if len==0
        throw(DivideError())
    end
    return abs(dot(plane.n, point .- plane.r)) / len
end

"""
Returns two 3-element vectors with the two possible probe-centers, if the probe touches the three atoms. 
Otherwise, (nothing, nothing) is returned. 
"""
function centerOfProbe(aw1::AtomWrapper{T}, aw2::AtomWrapper{T}, aw3::AtomWrapper{T}, probe_positions, probe_radius::T) where T

    if haskey(probe_positions, (aw1.index_in_list, aw2.index_in_list, aw3.index_in_list))
        value = probe_positions[aw1.index_in_list, aw2.index_in_list, aw3.index_in_list]
        if value!==nothing
            c1 = value.points[:, 1]
            c2 = value.points[:, 2]
            return c1, c2
        else
            return nothing, nothing
        end
    else
        s1 = deepcopy(aw1.atom)
        s2 = deepcopy(aw2.atom)
        s3 = deepcopy(aw3.atom)

        s1.radius += probe_radius
        s2.radius += probe_radius
        s3.radius += probe_radius

        point1, point2 = getIntersectionPoints(s1, s2, s3, false)

        if point1!==nothing && point2!==nothing
            pp = ProbePosition([ProbeStatus.NOT_TESTED, ProbeStatus.NOT_TESTED], [c1, c2])
            probe_positions[(aw1.index_in_list, aw2.index_in_list, aw3.index_in_list)] = pp
            return c1, c2
        else
            probe_positions[(aw1.index_in_list, aw2.index_in_list, aw3.index_in_list)] = nothing
            return nothing, nothing
        end
    end
end

function checkProbe(probe::ProbeSphere{T}, aw1, aw2, aw3, neighbourDict, all_atoms) where T
    idx1, idx2, idx3 = sort([aw1.index_in_list, aw2.index_in_list, aw3.index_in_list])
    pp = probe_positions[idx1, idx2, idx3]
    index = (probe.r == pp.points[:, 1]) ? 1 : 2

    if pp.status[index] == ProbeStatus.NOT_TESTED
        ok = true
        idx_list = neighboursOfThreeAtoms(aw1, aw2, aw3, neighbourDict)
        i = 1
        while ok && i<=length(idx_list)
            atom = all_atoms[idx_list].atom
            dist = (probe.radius + atom.radius)^2
            if isless_tolerance(norm(probe.r .- atom.r), dist)
                pp.status[index] = ProbeStatus.NOT_OK
                ok = false
            end
            i += 1
        end
        if ok
            pp.status[index] = ProbeStatus.OK
        end

    end
    return (pp.status[index] == ProbeStatus.OK)
end


"""
Wrapper function around the dictionary access (returns a Set of indices). 
It fills the dictionary the first time that a key is requested. 
"""
function neighboursOfTwoAtoms(aw1::AtomWrapper{T}, aw2::AtomWrapper{T}, neighbourDict::Dict{Tuple{Index, Index}, Set{Index}})

    if aw1.index_in_list>aw2.index_in_list
        temp = aw1
        aw1 = aw2
        aw2 = temp
    end

    if !haskey(neighbourDict, (aw1.index_in_list, aw2.index_in_list))
        shared_neighbours = intersect(aw1.neighbours, aw2.neighbours)
        neighbourDict[(aw1.index_in_list, aw2.index_in_list)] = shared_neighbours
        return shared_neighbours
    end

    return neighbourDict[(aw1.index_in_list, aw2.index_in_list)]
end

"""
Returns a Set of indices. 
"""
function neighboursOfThreeAtoms(aw1::AtomWrapper{T}, aw2::AtomWrapper{T}, aw3::AtomWrapper{T}, neighbourDict) where T
    set1 = neighboursOfTwoAtoms(aw1, aw2, neighbourDict)
    set2 = neighboursOfTwoAtoms(aw1, aw3, neighbourDict)

    return intersect(set1, set2)
end

"""
Computes probe positions for all atom_indices in neighbours_inds 
and returns a list of (atom_index, probe_sphere). 
"""
function findThirdAtom(aw1::AtomWrapper{T}, aw2::AtomWrapper{T}, neighbours::AbstractVect{AtomWrapper{T}}, probe_positions, probe_radius::T) where T
    # This function computes a list of all atoms (with its probe positions)
    # which can be touched by the probe sphere when it touches the two given
    # atoms 
    result = []
    for aw3 in neighbours
        c1, c2 = centerOfProbe(aw1, aw2, aw3, probe_positions, probe_radius)
        if c1!==nothing && c2!==nothing
            if !any(isnan.(c1))
                push!(result, (aw3.index_in_list, ProbeSphere(c1, probe_radius)))
            end
            if !any(isnan.(c2))
                push!(result, (aw3.index_in_list, ProbeSphere(c2, probe_radius)))
            end
        end
    end
    return result
end

function connectRSElements!(v1::RSVertex, v2::RSVertex, v3::RSVertex, e1::RSEdge, e2::RSEdge, e3::RSEdge, f, probe::ProbeSphere, v1Atom::Atom{T}, v2Atom::Atom{T}, v3Atom::Atom{T}) where T
    e1.vertices = (v1, v2)
    e1.faces[1] = f

    e2.vertices = (v2, v3)
    e2.faces[1] = f

    e3.vertices = (v3, v1)
    e3.faces[1] = f

    f.vertices = (v1, v2, v3)
    f.edges = (e1, e2, e3)

    plane = Plane(v1Atom.r, v2Atom.r, v3Atom.r)
    f.center = probe.r
    f.normal = plane.n
    if isless_tolerance(dot(f.normal, probe.r), dot(f.normal, v1Atom.r))
        f.normal .*= -1
    end
    f.singular = isless_tolerance(getDistance(probe.r, plane), probe.radius)
end

function insert(f::RSFace, rs::ReducedSurface, new_faces::Set{RSFace})
    insert(f, rs)
    push!(new_faces, f)

    push!(f.vertices[1].faces, f)
    push!(f.vertices[2].faces, f)
    push!(f.vertices[3].faces, f)

    for e in f.edges
        push!(e.vertices[1].edges, e)
        push!(e.vertices[2].edges, e)
    end
end

function insert(f::RSFace, rs::ReducedSurface)
    f.index = rs.number_of_faces
    rs.number_of_faces += 1
    push!(rs.faces, f)
end

function insert(v::RSVertex, rs::ReducedSurface, new_vertices::Set{RSVertex}, vertices)
    insert(v, rs)
    push!(new_vertices, v)
    push!(vertices[v.atom_index], v)
end

function insert(v::RSVertex, rs::ReducedSurface)
    v.index = rs.number_of_vertices
    rs.number_of_vertices += 1
    push!(rs.vertices, v)
end

function insert(e::RSEdge, rs::ReducedSurface, only_rs::Bool = false)
    e.index = rs.number_of_edges
    rs.number_of_edges += 1
    push!(rs.edges, e)

    if !only_rs
        push!(e.vertices[1].edges, e)
        push!(e.vertices[2].edges, e)
    end
end


function getRSComponent(faces::Vector{RSFace}, number_of_faces, new_faces::Set{RSFace}, all_atoms, neighbourdict, probe_radius, vertices, rs, rm_vertices, new_vertices, probe_positions)
    i=1
    while i<=number_of_faces
        if faces[i]!==nothing
            if !treatFace(faces[i], new_faces::Set{RSFace}, all_atoms, neighbourdict, probe_radius, vertices, rs, rm_vertices, new_vertices, probe_positions)
                i=1
            else
                i+=1
            end
        else
            i+=1
        end
    end
    extendComponent(new_vertices, rm_vertices, all_atoms, neighbourdict, probe_radius, new_faces, vertices, faces, number_of_faces, rs, probe_positions)
end

function treatFace(f::RSFace, new_faces::Set{RSFace}, all_atoms, neighbordict, probe_radius, vertices, rs, rm_vertices, new_vertices, probe_positions)
    if f.edges[1].faces[2]===nothing
        if !treatEdge(f.edges[1], all_atoms, neighbordict, probe_radius, vertices, rs, rm_vertices, new_vertices, new_faces, probe_positions)
            return false
        end
    end

    if f.edges[2].faces[2]===nothing
        if !treatEdge(f.edges[2], all_atoms, neighbordict, probe_radius, vertices, rs, rm_vertices, new_vertices, new_faces, probe_positions)
            return false
        end
    end

    if f.edges[3].faces[2]===nothing
        if !treatEdge(f.edges[3], all_atoms, neighbordict, probe_radius, vertices, rs, rm_vertices, new_vertices, new_faces, probe_positions)
            return false
        end
    end

    delete!(new_faces, f)
    return true
end

function treatEdge(edge::RSEdge, all_atoms, neighbordict, probe_radius, vertices, rs, rm_vertices, new_vertices, new_faces, probe_positions)
    # This function rolls the probe sphere over a RSEdge.
    # From all atoms that can be touched by the probe sphere when it touches
    # the two atoms of the edge is this one selected for which the rotation
    # angle is the smallest. A new face is found.
    # If this face already exists the edge exists twice, too. These two
    # edges and their vertices are joined.
    # If the face does not exist yet, it will be created. A new vertex and
    # two new edges will be created, too.
    # In both cases the treated edge will be updated. It has not to be
    # considerd again.

    # find third atom
    start_face = RSFace(edge->face_[0])		# the edge already knows the
    vertex1 = RSVertex(edge->vertex_[1])		# starting face and their
    vertex2 = RSVertex(edge->vertex_[2])		# two vertices
    idx1 = vertex1.atom_index
    idx2 = vertex2.atom_index
    try
        idx3, probe, phi = thirdAtom(vertex1,vertex2,start_face, all_atoms, neighbordict, probe_radius, vertices, rs, rm_vertices, probe_positions)
    catch e
        test_message = "PROBE SPHERE TOUCHES FOUR ATOMS"
        if e isa ErrorException && e.msg == test_message
            return false
        else
            rethrow(e)
        end
    end
    sphere1 = all_atoms[idx1].atom
    sphere2 = all_atoms[idx2].atom
    sphere3 = all_atoms[idx3].atom
    # build a new face and two new edges
    vertex3 = RSVertex(idx3)
    new_face = RSFace(  vertex1,vertex2,vertex3, 
                        nothing,nothing,nothing, 
                        probe.r,getFaceNormal(sphere1,sphere2,sphere3,probe),false,-1) # provisorial new face

    test = faceExists(new_face,vertices[vertex1.atom_index]) 
    if (test === nothing)
        # built face doesn't exist yet
        # The new vertex has to be created since we don't know at this time
        # whether it is a new vertex or not.
        # Attention: one atom can build more than one vertex!
        insert(vertex3, rs, new_vertices, vertices)
        all_atoms[idx3].status = AtomStatus.ON_SURFACE
        edge1 = RSEdge()
        edge1.vertices = (vertex2, vertex3)
        edge1.faces[1] = new_face

        edge2 = RSEdge()
        edge2.vertices = (vertex3, vertex1)
        edge2.faces[1] = new_face

        new_face.edges = (edge, edge1, edge2)
        TPlane3<double> plane(sphere1.r,sphere2.r,sphere3.r)
        new_face.singular = isless_tolerance(getDistance(probe.r,plane), probe_radius)
        insert(new_face, rs, new_faces)
    else
        # built face exitsts already
        # the corresponding edge in the existing face has to be found
        i, test_edge = getSimilarEdge(test, edge)
        if test_edge === nothing
            test_edge = RSEdge()
        end
        # Now the corresponding vertices of the corresponding edges have to be
        # joined and one of them has to be deleted (if they are not equal). This
        # is neccessary since creating a new face always creates a new vertex.
        test_vertex1 = test_edge.vertices[1]
        test_vertex2 = test_edge.vertices[2]
        if (test_vertex1.atom_index == vertex2.atom_index)
            tmp = test_vertex1
            test_vertex1 = test_vertex2
            test_vertex2 = tmp
        end
        # now we know which vertices are corresponding
        if vertex1 != test_vertex1
            # the vertices only have to be joined if they are not equal
            join!(vertex1, test_vertex1)
            substitute!(test_vertex1, vertex1)
            rs.vertices.[test_vertex1.index] =nothing
            push!(rm_vertices, test_vertex1)
            delete!(vertices.[test_vertex1.atom_index], test_vertex1)
        end
        if vertex2 != test_vertex2
            # the vertices only have to be joined if they are not equal
            join!(vertex2, test_vertex2)
            substitute!(test_vertex2, vertex2)
            rs.vertices[test_vertex2.index] = nothing
            push!(rm_vertices, test_vertex2)
            delete!(vertices[test_vertex2.atom_index], test_vertex2)
        end
        # The vertices should have only one of the two corresponding edges.
        # The other will be deleted later.
        delete!(vertex1.edges, test_edge)
        delete!(vertex2.edges, test_edge)
        # The face should have only one of the two corresponding edges, too.

        if i<i || i>3
            throw(BoundsError())
        end
        test.edges[i] = edge 
        # Now can we delete the build face and vertex and the double edge.
        if (test_edge.index_ != -1)		# this can happens after a correct step
            rs.edges[test_edge.index] = nothing
        end
        new_face = test
    end			# face exitsts test
    # update edge
    circle1,circle2,circle3 = getCircles(idx1,idx2, all_atoms)

    line = Line(sphere1.r, sphere2.r-sphere1.r)

    TVector3<double> ip1;		# intersection points between
    TVector3<double> ip2;		# the edge and the probe sphere
    ip1, ip2 = getIntersectionPoints(probe,line)
    if (ip1!==nothing && ip2!==nothing 
        && isless_tolerance(norm(ip1, sphere2.r), norm(ip2, sphere2.r))) # ip1 is the intersection point next to the firstvertex of the edge
        temp = ip1
        ip1 = ip2
        ip2 = temp								
    end
    edge.faces[2] = new_face
    edge.center_of_torus = circle1.r
    edge.radius_of_torus = circle1.radius
    edge.angle = phi
    edge.circle0 = circle2
    edge.circle1 = circle3
    edge.intersection_point0 = ip1
    edge.intersection_point1 = ip2
    edge.singular = ip1!==nothing && ip2!==nothing
    if (edge.index == -1)
        insert(edge, rs, true)
    end
    return true
end

function getFaceNormal(a1::Atom{T}, a2::Atom{T}, a3::Atom{T}, probe::ProbeSphere)
    plane = Plane(a1.r, a2.r, a3.r)
    normal = plane.n
    if isless_tolerance(dot(normal, probe.r), dot(normal, a1.r))
        normal .*= -1
    end
    return normal
end

"""
Returns (idx, probe, phi)
"""
function thirdAtom(vertex1, vertex2, face, all_atoms, neighbordict, probe_radius, vertices, rs, rm_vertices, probe_positions)
    # This function chooses from all atoms which can be touced by the probe
    # sphere when it touches the given two vertices this one, for which is
    # the rotation angle the smalest.
    # If the rotation angle equals zero, the probe sphere can touch four or
    # more atoms an an exception is thrown.
    # If no atom can be found an exception is thrown.
    idx1 = vertex1.atom_index
    idx2 = vertex2.atom_index
    idx_list = neighboursOfTwoAtoms(all_atoms[idx1], all_atoms[idx2], neighbordict)
    candidates = findThirdAtom(all_atoms[idx1], all_atoms[idx2], map(idx -> all_atoms[idx], idx_list), probe_positions, probe_radius)
    old_angle = 3*π
    two_pi = 2*π
    axis = all_atoms[idx1].atom.r .- all_atoms[idx2].atom.r
    test_vector = cross(face.normal, axis)
    idx3 = third(face, vertex1, vertex2).atom_index
    if isless_tolerance(dot(test_vector, all_atoms[idx3].atom.r), dot(test_vector, all_atoms[idx1].atom.r))
        axis .*= -1
    end
    sphere1 = deepcopy(all_atoms[idx1].atom)
    sphere2 = deepcopy(all_atoms[idx2].atom)
    sphere1.radius += probe_radius
    sphere2.radius += probe_radius
    circle = getIntersectionCircle(sphere1,sphere2)
    start_probe = face.center
    v1 = start_probe .- circle.r
    face_normal = face.normal
    
    third = []
    for k in candidates
        if (k[1] != idx3) || (k[2].r != start_probe) # not found the starting face
            v2 = k[2].r .- circle.r
            new_angle = getOrientedAngle(v1,v2,axis)
            if iszero_tolerance(new_angle) || isequal_tolerance(new_angle, two_pi)
                correct(k[1], vertices, rs, rm_vertices, all_atoms, probe_radius)
                throw(ErrorException("PROBE SPHERE TOUCHES FOUR ATOMS"))
            end
            if islessorqual_tolerance(new_angle, old_angle)
                if isless_tolerance(new_angle, old_angle)
                    old_angle = new_angle
                    for t in third
                        if (all_atoms[t[1]].status ==AtomStatus.UNKNOWN)
                            all_atoms[t[1]].status == AtomStatus.INSIDE
                        end
                    end
                    third = []
                end
                push!(third, k)
            else
                if (all_atoms[k[1]].status == AtomStatus.UNKNOWN)
                    all_atoms[k[1]].status == AtomStatus.INSIDE
                end
            end
        end
    end
    if length(third) > 1
        for k in third
            correct(k[1], vertices, rs, rm_vertices, all_atoms, probe_radius)
        end
        throw(ErrorException("PROBE SPHERE TOUCHES FOUR ATOMS"))
    end

    all_atoms[third[1][1]].status = AtomStatus.ON_SURFACE
    return third[1][1], third[1][2], old_angle
end

" Returns the third vertex of f that is neither v1 nor v2."
function third(face::RSFace, v1::RSVertex, v2::RSVertex)
    if f.vertices[1]==v1 || f.vertices[1]==v2
        if f.vertices[2]==v1 || f.vertices[2]==v2
            return f.vertices[3]
        else
            return f.vertices[2]
        end
    else
        return f.vertices[1]
    end
end

function getOrientedAngle(a, b, n)
    # Calculate the length of the two normals
    bl = sqrt(dot(a, a))
    el = sqrt(dot(b, b))
    bel = sqrt(dot(a, b))

    # if one or both planes are degenerated
    if (bl * el == 0)
      throw(DivideError())
	end
    bel /= (bl * el)
    if (bel > 1.0)
      bel = 1
    else if (bel < -1.0)
      bel = -1
    end

    acosbel = acos(bel)	# >= 0

    if ((          n[1] * (a[3] * b[2] - a[2] * b[3])
				 + n[2] * (a[1] * b[3] - a[3] * b[1])
				 + n[3] * (a[2] * b[1] - a[1] * b[2])) > 0)
    	acosbel = 2*π - acosbel
	end

    return acosbel
end

function correct(idx, vertices::Vector{Vector{RSVertex}}, rs::ReducedSurface, rm_vertices, all_atoms, probe_radius)
    treat_faces = Set{RSFace}()
    test_vertices = Set{RSVertex}()
    delete_edges = Set{RSEdge}()

    for vertex in vertices[idx]
        empty!(treat_faces)
        empty!(test_vertices)
        empty!(delete_edges)

        faces = vertex.faces
        for f in faces
            remove!(f, delete_edges,test_vertices,treat_faces)
        end

        for f in faces
            delete!(treat_faces, f)
            delete!(new_faces, f)
            rs.faces[f.index] = nothing
        end

        for f in treat_faces
            rs.faces[f.index] = nothing
            push!(rs.faces, f)
            f.index = rs.number_of_faces
            rs.number_of_faces += 1
        end

        for edge in delete_edges
            index = edge.index
            if (index != -1)
                rs.edges[index] = nothing
            end
        end

        delete!(test_vertices, vertex)
        for test in test_vertices
            if length(test.edges)==0
                rs.vertices[test.index] = nothing
                delete!(vertices[test.atom_index], test)
                push!(rm_vertices, test)
            end
        end
        rs.vertices[vertex.index] = nothing
        delete!(vertices[idx], vertex)
        push!(rm_vertices, vertex)
    end
    rs.atoms[idx].atom.radius -= 0.001
    rs.atoms[idx].status = AtomStatus.UNKNOWN
    correctProbePosition(idx, probe_radius, all_atoms)
end

function remove!(face, del_edges, del_vertices, del_faces)
    delete!(face.vertices[1].faces, f)
    delete!(face.vertices[2].faces, f)
    delete!(face.vertices[3].faces, f)
    if face.edges[1] !== nothing
        if face.edges[1].faces[2] === nothing
            delete!(face.edges[1].vertices[1].edges, face.edges[1])
            delete!(face.edges[1].vertices[2].edges, face.edges[1])
            push!(del_vertices, faces.edges[1].vertices[1])
            push!(del_vertices, faces.edges[1].vertices[2])
            push!(del_edges, faces.edges[1])
            face.edges[1] = nothing
        else
            if(face.edges[1].faces[2] == f)
                face.edges[1].faces[2] = nothing
            elseif face.edges[1].faces[1] == f
                face.edges[1].faces[1] = face.edges[1].faces[2]
                face.edges[1].faces[2] = nothing
            end
            push!(del_faces, f)
        end
    end
    if face.edges[2] !== nothing
        if face.edges[2].faces[2] === nothing
            delete!(face.edges[2].vertices[1].edges, face.edges[2])
            delete!(face.edges[2].vertices[2].edges, face.edges[2])
            push!(del_vertices, face.edges[2].vertice[1])
            push!(del_vertices, face.edges[2].vertice[2])
            push!(del_edges, face.edges[2])
            face.edges[2] = nothing
        else
            if(face.edges[2].faces[2] == f)
                face.edges[2].faces[2] = nothing
            elseif face.edges[2].faces[1] == f
                face.edges[2].faces[1] = face.edges[2].faces[2]
                face.edges[2].faces[2] = nothing
            end
            push!(del_faces, f)
        end
    end
    if face.edges[3] !== nothing
        if (face.edges[3].faces[1] === nothing)
            delete!(face.edges[3].vertices[1].edges, face.edges[3])
            delete!(face.edges[3].vertices[2].edges, face.edges[3])
            push!(del_vertices, face.edges[3].vertice[1])
            push!(del_vertices, face.edges[3].vertice[2])
            push!(del_edges, faces.edges[3])
            faces.edges[3] = nothing
        else
            if(face.edges[3].faces[2] == f)
                face.edges[3].faces[2] = nothing
            elseif face.edges[3].faces[1] == f
                face.edges[3].faces[1] = face.edges[3].faces[2]
                face.edges[3].faces[2] = nothing
            end
            push!(del_faces, f)
        end
    end
end


function correctProbePosition(idx::Index, probe_positions, probe_radius, all_atoms)
    for (key, value) in probe_positions
        if ((key[1] == idx) || (key[2] == idx) || (key[3] == idx))
            correctProbePosition(key, probe_positions, probe_radius, all_atoms)
        end
    end
end

function correctProbePosition(idx::NTuple{3, Index}, probe_positions, probe_radius, all_atoms)

    s1 = deepcopy(all_atoms[idx[1]])
    s1.atom.radius += probe_radius    
    s2 = deepcopy(all_atoms[idx[2]])
    s2.atom.radius += probe_radius
    s3 = deepcopy(all_atoms[idx[3]])
    s3.atom.radius += probe_radius

    c1, c2 = getIntersectionPoints(s1, s2, s3)
    if c1!==nothing && c2!==nothing
        if !haskey(probe_positions, idx)
            position = ProbePosition()
            probe_positions[idx] = position
        else
            position = probe_positions[idx]
        end
        position.status = [ProbeStatus.NOT_TESTED, ProbeStatus.NOT_TESTED]
        position.points[:, 1] = c1
        position.points[:, 2] = c2
    else
        delete!(probe_positions, idx)
    end
end

function faceExiste(face, vertices)
    for v in vertices
        if face ∈ v.faces
            return face
        end
    end
    return nothing
end

"Returns (index, similarEdge) or (-1, nothing)"
function getSimilarEdge(face::RSFace, edge::RSEdge)
    for i=1:3
        if isSimilar(face.edges[i], edge)
            return i, face.edges[i]
        end
    end
    return -1, nothing
end

function isSimilar(e1::RSEdge, e2::RSEdge)
    return (
        ((e1.vertices[1].atom_index == e2.vertices[1].atom_index) &&
         (e1.vertices[2].atom_index == e2.vertices[2].atom_index)    ) ||
        ((e1.vertices[1].atom_index == e2.vertices[2].atom_index) &&
         (e1.vertices[2].atom_index == e2.vertices[1].atom_index)    )
     )
end

function isSimilar(v1::RSVertex, v2::RSVertex)
    return v1.atom_index == v2.atom_index
end

"Adds edges and faces from v2 to v1. Returns a success-boolean. "
function join!(v1::RSVertex, v2::RSVertex)
    if isSimilar(v1, v2)
        union!(v1.edges, v2.edges)
        union!(v1.faces, v2.faces)
        return true
    end
    return false
end

"Substitues vertex references in edges and faces of v1 with v2. "
function substitute!(v1::RSVertex, v2::RSVertex)
    if isSimilar(v1, v2)
        for e in v1.edges
            if e.vertices[1]==v1
                e.vertices[1] = v2
            else if e.vertices[2]==v1
                e.vertices[2] = v2
            end
        end

        for f in v1.faces
            for i=1:3
                if f.vertices[i]==v1
                    f.vertices[i] = v2
                end
            end
        end
        return true
    end
    return false
end

"Return 3 circles or (nothing, nothing, nothing)"
function getCircles(idx1, idx2, all_atoms::AbstractVector{Atom{T}}) where T
    a1 = deepcopy(all_atoms[idx1].atom)
    a2 = deepcopy(all_atoms[idx2].atom)

    a1.radius += probe_radius
    a2.radius += probe_radius

    c1 = getIntersectionCircle(a1, a2)
    if c1!==nothing
        ratio = all_atoms[idx1].atom.radius / a1.radius

        c2 = Circle(a1.r + (c1.r-a1.r)*ratio, Vector{T}(3, undef), c1.radius * ratio)

        ratio = all_atoms[idx2].atom.radius / a2.radius
        c3 = Circle(a2.r + (c1.r-a2.r)*ratio, Vector{T}(3, undef), c1.radius * ratio)

        return c1, c2, c3
    end
    return nothing, nothing, nothing
end

function extendComponent(new_vertices, rm_vertices, all_atoms, neighbourdict, probe_radius, new_faces, vertices, faces, number_of_faces, rs, probe_positions)
    while !empty(new_vertices)
        face = nothing
        vertex1 = popfirst!(new_vertices)

        # Caution: vertex1 might be a dangling pointer!
        # Skip iteration if vertex has been deleted.
        if vertex1 ∈ rm_vertices
            continue
        end

        idx1 = vertex1.atom_index
        stop = false
        for i in all_atoms[idx1].neighboursOfThreeAtoms
            if stop
                break
            end
            if (all_atoms[i].status == AtomStatus.UNKNOWN)
                idx2 = i
                idx_list = neighboursOfTwoAtoms(all_atoms[idx1], all_atoms[idx2], neighbourdict)
                candidates = findThirdAtom(all_atoms[idx1], all_atoms[idx2], map(idx -> all_atoms[idx], idx_list), probe_positions, probe_radius)
                if (empty(candidates))
                    vertex2 = RSVertex(idx2)
                    edge = createFreeEdge(vertex1,vertex2, all_atoms, probe_radius, neighbourdict)
                    if (edge !== nothing)
                        insert(edge, rs)
                        insert(vertex2, rs, new_faces, vertices)
                        all_atoms[idx2].status = AtomStatus.ON_SURFACE
                        push!(new_vertices, vertex1)
                        push!(new_vertices, vertex2)
                        # i = neighbours_[atom1].end()--; ???
                        break
                    end
                else
                    for (idxJ, probe) in candidates
                        if (all_atoms[idxJ].status == AtomStatus.UNKNOWN)
                            if (checkProbe(probe, all_atoms[idx1], all_atoms[idx2], all_atoms[idxJ], neighbourdict, all_atoms))
                            
                                face = RSFace()
                                edge1 = RSEdge()
                                edge2 = RSEdge()
                                edge3 = RSEdge()
                                vertex2 = RSVertex(idx2)
                                vertex3 = RSVertex(idxJ)
                                connectRSElements!(vertex1,vertex2,vertex3,
                                                                     edge1,edge2,edge3,
                                                                     face,probe, all_atoms[idx1], all_atoms[idx2], all_atoms[idxJ]);
                                insert(face, rs, new_faces)
                                insert(vertex2, rs, new_vertices, vertices)
                                insert(vertex3, rs, new_vertices, vertices)
                                all_atoms[idx2].status = AtomStatus.ON_SURFACE
                                all_atoms[idxJ].status = AtomStatus.ON_SURFACE
                                push!(new_vertices, vertex1)
                                push!(new_vertices, vertex2)
                                push!(new_vertices, vertex3)
                                #i = neighbours_[atom1].end()--;
                                #j = candidates.end()--;
                                #????
                                stop = true
                                break
                            end
                        end
                    end
                end
            end
        end
        if face !== nothing
            getRSComponent(faces, number_of_faces, new_faces, all_atoms, neighbordict, probe_radius, vertices, rs, rm_vertices, new_vertices, probe_positions)
        end
    end
    empty!(rm_vertices)
end

function createFreeEdge(v1::RSVertex, v2::RSVertex, all_atoms, probe_radius, neighbourdict)
    idx1 = v1.atom_index
    idx2 = v2.atom_index

    # compute the three circles describing the toric face
    circle1, circle2, circle3 = getCircles(idx1, idx2, all_atoms)
    if ( !any(isnan.([circle1, circle2, circle3])) && # the probe hulls intersect
            isgreater_tolerance(circle1.radius, probe_radius))   # the radius of the toric edge is > 0 
        plane = Plane(circle1.r, circle1.n)

        TCircle3<double> test_circle;

        idx_list = neighboursOfTwoAtoms(all_atoms[idx1], all_atoms[idx2], neighbourdict);

        # find the mutual neighbours of both atoms
        for i in idx_list
            # put a sphere into the neighbour
            sphere = ProbeSphere(all_atoms[i].atom.r, all_atoms[i].atom.radius+probe_radius);

            test_circle = getIntersectionCircle(sphere,plane)
            if test_circle!==nothing
                radius_dist = test_circle.radius-circle1.radius
                radius_sum  = test_circle.radius+circle1.radius
                center_dist = norm(test_circle.r - circle1.r) # TODO find getSquareDistance

                if (islessorqual_tolerance(radius_dist*radius_dist, center_dist) 
                        && isgreaterorequal_tolerance(radius_sum*radius_sum, center_dist) ) # the circles intersect
                    return nothing
                end
            end
        end
        v::Vector{T} = [0, 0, 0]
        edge = RSEdge((vertex1, vertex2), (nothing, nothing), -1, circle1.r, circle1.radius,
                                    2*π, circle2, circle3,
                                    v, v, false)

        return edge
    end

    return nothing
end

"""
Solves x^2 + bx + c = 0
Returns (num_solutions, x1, x2). 
When num_solutions==0, x1 and x2 are nothing. When num_solutions==1, x1 is duplicated in x2. 
"""
function SolveQuadraticEquation(a, b, c)
    if a==0
        if b==0
            return 0, nothing, nothing
        end 
        return 1, c/b, c/b
    end

    discriminant = b^2 -4*a*c
    if isless_tolerance(discriminant, 0)
        return 0, nothing, nothing
    end

    sqrt_discriminant = sqrt(discriminant)
    if iszero_tolerance(sqrt_discriminant)
        return 1, -b/(2*a), -b/(2*a)
    else
        return 2, -b + sqrt_discriminant/(2*a), -b - sqrt_discriminant/(2*a)
    end
end
        