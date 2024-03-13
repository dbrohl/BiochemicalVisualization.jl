const Index = Int

@enumx AtomStatus begin
    ON_SURFACE,
    INSIDE,
    UNKNOWN
end

@enumx ProbeStatus begin
    OK,
    NOT_OK,
    NOT_TESTED
end

struct AtomWrapper{T}
    atom::Atom{T}

    status::AtomStatus
    neighbours::Set{Index}

    index_in_list::Index

    function AtomWrapper(atom::Atom{T}) where T
        new{T}(atom, AtomStatus.UNKNOWN, Set(), 0)
    end

end

struct Circle{T}
    r::Vector{T} # center
    n::Vector{T} # normal
    radius::T
end

struct Line{T}
    r::Vector{T}
    d::Vector{T}

    Line{T}() = new{T}(Vector{T}(undef, 3), Vector{T}(undef, 3))
end

struct Plane{T}
    r::Vector{T} # point 
    n::Vector{T} # normal

    function Plane{T}(n::Vector{T}, d::T) where T
        r = Vector{T}(undef, 3)
        if all(iszero_tolerance.(n))
            throw(DivideError())
        end
        if !iszero_tolerance(n[1])
            r = [-d / n[1], T(0), T(0)]
        elseif !iszero_tolerance(n[2])
            r = [T(0), -d / n[2], T(0)]
        elseif !iszero_tolerance(n[3])
            r = [T(0), T(0), -d / n[3]]
        end
        return Plane{T}(r, n)
    end

    function Plane{T}(a::vector{T}, b::vector{T}, c::vector{T})
        return Plane{T}(a, cross(a.-b, b.-c))
    end
end

mutable struct RSVertex
    atom_index::Index
    edges::Set{RSEdge}
    faces::Set{RSFace}

    index::Index

    RSVertex(index) = new(index, Set(), Set())

    RSVertex(v::RSVertex) = new(v.atom_index, v.edges, v.faces, v.index)
end

mutable struct RSEdge{T}
    vertices::NTuple{2, Union{Nothing, RSVertex}}
    faces::NTuple{2, Union{Nothing, RSFace}}

    index = -1

    center_of_torus::Union{Nothing, Vector{T}}
    radius_of_torus::Union{Nothing, T}
    angle::Union{Nothing, T}
    circle0::Union{Nothing, Circle{T}}
    circle1::Union{Nothing, Circle{T}}
    intersection_point0::Union{Nothing, Vector{T}}
    intersection_point1::Union{Nothing, Vector{T}}
    singular::Union{Nothing, Bool}
    RSEdge{T}() = new{T}((nothing, nothing), (nothing, nothing), -1,
    nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end

mutable struct RSFace{T}
    vertices::NTuple{3, Union{Nothing, RSVertex}}
    edges::NTuple{3, Union{Nothing, RSEdge}}
    center::Vector{T}
    normal::Vector{T}
    singular::Bool

    index::Index

    function RSFace{T}()
        part = new{T}()
        part.vertices = (nothing, nothing, nothing)
        part.edges = (nothing, nothing, nothing)
        part.index = -1
        return part
    end

    function RSFace{T}(f::RSFace{T}) where T
        new{T}(f.vertices, f.edges, f.center, f.normal, f.singular, f.index)
    end

    function RSFace{T}(v1, v2, v3, e1, e2, e3, center, normal, singular, index)
        new{T}((v1, v2, v3), (e1, e2, e3), center, normal, singular, index)
end

struct ProbeSphere{T}
    r::Vector{T}
    radius::T
end

struct ProbePosition{T}
    status::Vector{ProbeStatus.T}
    points::Matrix{T} # (3, 2) matrix
end


struct RSComputerStore
    # Notes: Dicts can be handled as arrays with key==index. The key is the index of the atom in the large atom-array

    # Properties that are related to a single atom can be stored directly in the atom
    atom_status::Dict{Index, AtomStatus.T}
    neighbours::Dict{Index, Set{Index}} # for each atom a list of its neighbours

    probe_radius::T
    neighboursOfTwoAtoms::Dict{Tuple{Index, Index}, Set{Index}} # for each pair of atoms (indices are ordered so that ind_a<ind_b) a list of their common neighbours
    probePositions::Dict{Tuple{Index, Index, Index}, Union{Nothing, ProbePosition}}

    newVertices::Set{RSVertex}
    newFaces::Set{RSFace}

    vertices::Vector{Vector{RSVertex}}

    rm_vertices::Set{RSVertex}
end

struct ReducedSurface
    number_of_vertices::Index
    number_of_edges::Index
    number_of_faces::Index

    atoms::Vector{AtomWrapper{T}}

    vertices::Vector{RSVertex}
    edges::Vector{Union{Nothing, RSEdge}}
    faces::Vector{Union{Nothing, RSFace}} 
end