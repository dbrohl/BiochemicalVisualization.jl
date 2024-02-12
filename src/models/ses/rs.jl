# create reduced surface
    # create list of spheres(atom positions and radii+probe_radius (+ ggf. offset))
    # preprocessing
        # determine bounding box of atom centers (without vdw/probe radius)
        # determine max radius
        # create hash grid (Sanner et al use a tree instead of a grid. )

#TODO separate methode für chains? sollten chains als verschiedene Oberflächen behandelt werden?


function create_reduced_surface(
    ac::System{T}, probe_radius::T) where {T<:Real}
    if natoms(ac)==0
        return # TODO better error
    end

    # void RSComputer::preProcessing()
    epsilon = 1e-4;

    first_atom = first(eachatom(ac))
    max_radius = 0
    bounding_box = [first_atom.r;; first_atom.r] # [x/y/z, min/max]

    for atom in eachatom(ac)
        max_radius = max(max_radius, atom.radius==0 ? 1 : atom.radius) #TODO why are radii sometimes zero? how should this be handled in special cases?

        bounding_box[:, 1] = min.(bounding_box[:, 1], atom.r)
        bounding_box[:, 2] = max.(bounding_box[:, 2], atom.r)
    end
    log_info(rs, "Bounding box of all $(natoms(ac)) atoms", bounding_box)

    max_dist = 2*(max_radius + probe_radius)
    hashgrid = HashGrid{T, Atom{T}}(bounding_box, [max_dist, max_dist, max_dist], atom -> atom.r, natoms(ac)) #TODO eventuell bounding_box anpassen und vergrößern (->reducedSurface.C, line 1861f.)
    # TODO ggf. mit CellListMap machen, falls das schneller ist

    for atom in eachatom(ac)
        too_close = false
        atom_radius = atom.radius==0 ? 1 : atom.radius
        for (neighbor, _) in each_neighbor(hashgrid, atom.r)
            neighbor_radius = neighbor.radius==0 ? 1 : neighbor.radius

            if (norm(atom.r .- neighbor.r) - 0.05*max(neighbor_radius, atom_radius)) < epsilon

                too_close = true
                break
            end
        end

        if !too_close
            push!(hashgrid, atom)
        end
    end
    log_info(rs, "Hashgrid contains $(length(hashgrid.data)) atoms")


    # TODO cache actual neighbors (outside of the hashgrid since they are dependant on individual radii)
    # why sorting? can i store them inside of atom?
    # for atom in eachatom(ac) #TODO bis auf das letzte
    #     for (neighbor, inter_atom_dist) in each_neighbor(hashgrid, atom.r, max_dist)
    #         if neighbor.number > atom.number && 

    start = 1;
    while start!=0
        start = getStartPosition()

        if start==2
            extendComponent()
        elseif start==3
            getRSComponent()
        end
    end

    rs_->clean();

    #TODO if Constants::EPSILON is used later, it has to be the old value from line 676 in reducedSurface.C



end


function getStartPosition()
    if findFirstFace()!==nothing
        return 3
    end
    if findFirstEdge()!==nothing
        return 2
    end
    if findFirstVertex()!==nothing
        return 1
    end
    return 0
end

function findFirstFace()
    for direction = 1:3
        for extreme = 1:2
            face = findFace(direction, extreme)
            if face!==nothing
                return face
            end
        end
    end
    return nothing
end
function findFirstEdge()
    for direction = 1:3
        for extreme = 1:2
            edge = findEdge(direction, extreme)
            if edge!==nothing
                return edge
            end
        end
    end
    return nothing
end
function findFirstVertex(grid::HashGrid)
    for atom in grid.data
        if atom.status==UNKNOWN # TODO make status a property of the contianers in grid
            if length(atom.neighbors)==0
                # TODO
                # RSVertex* vertex = new RSVertex(i);
                # insert(vertex);
                
                # return vertex;
            end
        end
    end
    return nothing
end

function findFace(direction, extreme, grid::HashGrid{AtomWrapper})
    a1 = findFirstAtom(direction, extrem, grid.data)
    if a1 === nothing
        return nothing
    end

    a2 = findSecondAtom(a1, direction, extrem, grid)
    if a2 === nothing
        return nothing
    end

    const std::deque<Index>& s = neighboursOfTwoAtoms(SortedPosition2(a1, a2));
    candidates = findThirdAtom(a1, a2, s);
    if candidates.empty()
        return nothing
    end
    
    std::deque<std::pair<Index,TSphere3<double> > >::iterator i = candidates.begin();
    Index a3 = -1;
    TSphere3<double> probe;
    bool found = false;
    while (!found && i != candidates.end())
    {
        a3 = i->first;
        probe = i->second;
        found = (atom_status_[a3] == STATUS_UNKNOWN) &&
                        checkProbe(probe,SortedPosition3(a1,a2,a3));
        i++;
    }
    if (found)
    {
        RSVertex* vertex1 = new RSVertex(a1);
        RSVertex* vertex2 = new RSVertex(a2);
        RSVertex* vertex3 = new RSVertex(a3);
        
        RSEdge* e1 = new RSEdge;
        RSEdge* e2 = new RSEdge;
        RSEdge* e3 = new RSEdge;
        
        RSFace* face = new RSFace;
        
        updateFaceAndEdges(vertex1,vertex2,vertex3,e1,e2,e3,face,probe);
        
        insert(face);
        insert(vertex1);
        insert(vertex2);
        insert(vertex3);
        
        return face;
    }
    else
    {
        atom_status_[a1] = STATUS_INSIDE;
        atom_status_[a2] = STATUS_INSIDE;
        return NULL;
    }

end

function findFirstAtom(direction, extreme, atoms::AbstractVector{AtomWrapper})
    extreme_atom_index = nothing
    # find the first atom of unknown status
    index = nothing
    for (i, atom) in enumerate(atoms)
        if atom.status == UNKNOWN
            index = i
            break
        end
    end

    
    if index!==nothing
        extreme_atom_index = index
        next_atom = atoms[index]
        extreme_value = ((extreme == 1) ? next_atom.atom.r[direction]-next_atom.atom.radius 
                                        : next_atom.atom.r[direction]+next_atom.atom.radius)

        #find the atom of unknown status lying on the extrem position
        for i=index+1:length(atoms)
            if atoms[i].status == UNKNOWN
                next_atom = atoms[i]
                temp = ((extreme == 0) ? next_atom.atom.r[direction]-next_atom.atom.radius
                                       : next_atom.atom.r[direction]+next_atom.atom.radius)
                if (((extreme == 0) && isless_tolerance(temp, extreme_value)) ||
                        ((extreme != 0) && isgreater_tolerance(temp, extreme_value)))
                    extreme_value = temp
                    extreme_atom_index = i
                end
            end
        end
    end
    return extreme_atom_index
end

function findSecondAtom(atom_index_1, direction, extreme, grid, probe_radius)
    second_atom_index = nothing;
    # find the first neighbour atom of unknown status
    selected_atom = nothing
    iter = each_neighbor(grid, atoms[atom_index_1].atom.r)
    for (neighbor, distance) in iter
        if neighbor.status == UNKNOWN
            selected_atom = neighbor
            break
        end
    end



    if selected_atom!==nothing
        second_atom_index = findfirst(x->x==selected_atom, grid.data)
        first_atom = atoms[atom_index_1]
        first_atom.atom.radius += probe_radius;
        extreme_value = ((extreme == 0) ? first_atom.atom.r[direction]+first_atom.atom.radius
                                        : first_atom.atom.r[direction]-first_atom.atom.radius)
        next_atom = nothing
        # find the neighbour atom of unknown status lying on the extreme position
        for (neighbor, distance) in iter
            if neighbor.status == UNKNOWN
                next_atom = neighbor
                next_atom.radius += probe_radius; # ensure that this doesnt destroy the underlying data
                intersection_circle = getIntersection(first_atom,next_atom)
                if intersection_circle!==nothing
                    next_extreme = getCircleExtremum(intersection_circle,direction,extreme);
                    if (((extreme == 0) && isless_tolerance(next_extreme,extreme_value)) ||
                            ((extreme != 0) && isgreater_tolerance(next_extreme,extreme_value)))
                        extreme_value = next_extreme;
                        second_atom_index = findfirst(x->x==neighbor, grid.data);
                    end
                end
            end
        end
    end
    return second_atom_index
end

function findThirdAtom(atom_index_1, atom_index_2, third_candidates, all_atoms)
    # This function computes a list of all atoms (with its probe positions)
    # which can be touched by the probe sphere when it touches the two given
    # atoms
    std::pair<Index, TSphere3<double> > candidate;
    TVector3<double> center1, center2;
    TSphere3<double> probe;
    probe.radius = rs_->probe_radius_;

    result = []

    for (i, atom) in enumerate(third_candidates)
        if (centerOfProbe(SortedPosition3(atom1,atom2,*i),center1,center2))
            if (!(Maths::isNan(center1.x) || Maths::isNan(center1.y) || Maths::isNan(center1.z)))
                probe.p.set(center1);
                candidate.first = *i;
                candidate.second = probe;
                push!(result, candidate);
            end

            if (!(Maths::isNan(center2.x) || Maths::isNan(center2.y) || Maths::isNan(center2.z)))
                probe.p.set(center2);
                candidate.first = *i;
                candidate.second = probe;
                push!(result, candidate);
            end
        end
    end
    return result
end

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
    if isgreaterorequals_tolerance(abs(a.radius - b.radius), dist)
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

function getIntersectionPoints(s1::Atom{T}, s2::Atom{T}, s3::Atom{T}, test=true)

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
    TLine3<T> line;
    if (GetIntersection(plane1,plane2,line))
    {
        TVector3<T> diff(s1.p-line.p);
        T x1 = 0;
        T x2 = 0;
        if (SolveQuadraticEquation(line.d*line.d, -diff*line.d*2, diff*diff-r1_square, x1,x2) > 0)
        {
            p1 = line.p+x1*line.d;
            p2 = line.p+x2*line.d;
            if (test)
            {
                TVector3<T> test = s1.p-p1;
                if (Maths::isNotEqual(test*test,r1_square))
                {
                    return false;
                }
                test = s1.p-p2;
                if (Maths::isNotEqual(test*test,r1_square))
                {
                    return false;
                }
                test = s2.p-p1;
                if (Maths::isNotEqual(test*test,r2_square))
                {
                    return false;
                }
                test = s2.p-p2;
                if (Maths::isNotEqual(test*test,r2_square))
                {
                    return false;
                }
                test = s3.p-p1;
                if (Maths::isNotEqual(test*test,r3_square))
                {
                    return false;
                }
                test = s3.p-p2;
                if (Maths::isNotEqual(test*test,r3_square))
                {
                    return false;
                }
            }
            return true;
        }
    }
    return false;

    return p1, p2

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

function centerOfProbe(i_a, i_b, i_c)
    back = false
    if (i_a, i_b, i_c) ∈ probe_positions.keys()
        if probe_positions[i_a, i_b, i_c]!==nothing
            c1 = probe_positions[i_a, i_b, i_c].point[1]
            c2 = probe_positions[i_a, i_b, i_c].point[2]
            back = true
        end
    else
        s1 = atoms[i_a]
        s2 = atoms[i_a]
        s3 = atoms[i_a]

        s1.radius += probe_radius
        s2.radius += probe_radius
        s3.radius += probe_radius

        result1, result2 = getIntersectionPoints(s1, s2, s3, false)

        if ()
        {
            ProbePosition* position = new ProbePosition;
            position->status[0] = STATUS_NOT_TESTED;
            position->status[1] = STATUS_NOT_TESTED;
            position->point[0] = c1;
            position->point[1] = c2;
            probe_positions_.insert(std::make_pair(pos, position));
            back = true;
        }
        else
        {
            probe_positions_.insert(std::make_pair(pos, (ProbePosition*)NULL));
        }
        end

    return back;
end