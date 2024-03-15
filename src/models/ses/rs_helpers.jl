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

function map(f, s::Set, resultType)
    result = Set{resultType}()
    for value in s
        push!(result, f(value))
    end
    return result
end

# ----- geometry functions -----

""" Returns the intersection circle or nothing. """
function getIntersectionCircle(a::Atom{T}, b::Atom{T}) where T
    norm = b.r .- a.r
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

    return Circle{T}(a.r .+ (norm .* length), norm / dist, sqrt(square_radius))
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
function getIntersectionPoints(s1::Atom{T}, s2::Atom{T}, s3::Atom{T}, test::Bool=true) where T

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
        if e isa DivideError
            log_info(surface, "plane was not created", s2.r .- s1.r, s3.r .- s1.r)
            return nothing, nothing
        else 
            rethrow(e)
        end
    end

    line = getIntersectionLine(plane1,plane2)
    if line!==nothing
        diff = s1.r .- line.r
        num_solutions, x1, x2 = SolveQuadraticEquation(dot(line.d, line.d), dot(-diff, line.d)*2, dot(diff, diff)-r1_square)
        if (num_solutions > 0)
            p1 = line.r .+ x1.*line.d
            p2 = line.r .+ x2.*line.d
            if test
                testVec = s1.r .- p1
                if isnotequal_tolerance(test^2, r1_square)
                    return nothing, nothing
                end
                testVec = s1.r .- p2
                if isnotequal_tolerance(test^2, r1_square)
                    return nothing, nothing
                end
                testVec = s2.r .- p1
                if isnotequal_tolerance(test^2, r2_square)
                    return nothing, nothing
                end
                testVec = s2.r .- p2
                if isnotequal_tolerance(test^2, r2_square)
                    return nothing, nothing
                end
                testVec = s3.r .- p1
                if isnotequal_tolerance(test^2, r3_square)
                    return nothing, nothing
                end
                testVec = s3.r .- p2
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

    return line.r .+ x1 .* line.d, line.r .+ x2 .* line.d
end

"""
Computes the intersection line between plane1 and plane2. 
Returns either a Line or nothing. 
"""
function getIntersectionLine(plane1::Plane{T}, plane2::Plane{T}) where T
    u = dot(plane1.r, plane1.n)
    v = dot(plane2.r, plane2.n)
    det = plane1.n[1]*plane2.n[2] - plane1.n[2]*plane2.n[1]
    line = Line{T}()
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

function getCircleExtremum(circle::Circle, direction, extreme)
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
    dict_key = Tuple(sort([aw1.index_in_list, aw2.index_in_list, aw3.index_in_list]))

    if haskey(probe_positions, dict_key)
        value = probe_positions[dict_key]
        if value!==nothing
            c1 = value.points[1]
            c2 = value.points[2]
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
            pp = ProbePosition([ProbeStatus.NOT_TESTED, ProbeStatus.NOT_TESTED], [point1, point2])
            probe_positions[dict_key] = pp
            return point1, point2
        else
            probe_positions[dict_key] = nothing
            return nothing, nothing
        end
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
    elseif (bel < -1.0)
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

"Return 3 circles or (nothing, nothing, nothing)
circle2 and circle3 do not have a normal!
"
function getCircles(idx1, idx2, all_atoms::AbstractVector{AtomWrapper{T}}, probe_radius) where T
    a1 = deepcopy(all_atoms[idx1].atom)
    a2 = deepcopy(all_atoms[idx2].atom)

    a1.radius += probe_radius
    a2.radius += probe_radius

    c1 = getIntersectionCircle(a1, a2)
    if c1!==nothing
        ratio = all_atoms[idx1].atom.radius / a1.radius

        c2 = Circle(Vector{T}(a1.r + (c1.r-a1.r)*ratio), Vector{T}(undef, 3), c1.radius * ratio)

        ratio = all_atoms[idx2].atom.radius / a2.radius
        c3 = Circle(Vector{T}(a2.r + (c1.r-a2.r)*ratio), Vector{T}(undef, 3), c1.radius * ratio)

        return c1, c2, c3
    end
    return nothing, nothing, nothing
end

# ----- insert functions -----

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


# ----- convenience for RSX structs -----

"Returns (index, similarEdge) or (0, nothing)"
function getSimilarEdge(face::RSFace, edge::RSEdge)
    for i=1:3
        if isSimilar(face.edges[i], edge)
            return i, face.edges[i]
        end
    end
    return 0, nothing
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
            elseif e.vertices[2]==v1
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

function faceExists(face, vertices)
    for v in vertices
        if face ∈ v.faces
            return face
        end
    end
    return nothing
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

function getFaceNormal(a1::Atom{T}, a2::Atom{T}, a3::Atom{T}, probe::ProbeSphere) where T
    plane = Plane(a1.r, a2.r, a3.r)
    normal = plane.n
    if isless_tolerance(dot(normal, probe.r), dot(normal, a1.r))
        normal .*= -1
    end
    return normal
end

# ----- neighbour sets -----
"""
Wrapper function around the dictionary access (returns a Set of indices). 
It fills the dictionary the first time that a key is requested. 
"""
function neighboursOfTwoAtoms(aw1::AtomWrapper{T}, aw2::AtomWrapper{T}, neighbourDict::Dict{Tuple{Index, Index}, Set{Index}}) where T

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