# abstract types for recursive data structures
abstract type AbstractVertex end
abstract type AbstractEdge end
abstract type AbstractFace end

const Index = Int

@enumx AtomStatus begin
    ON_SURFACE
    INSIDE
    UNKNOWN
end

@enumx ProbeStatus begin
    OK
    NOT_OK
    NOT_TESTED
end

mutable struct AtomWrapper{T}
    atom::Atom{T}

    status::AtomStatus.T
    neighbours::Set{Index}

    index_in_list::Index

    function AtomWrapper(atom::Atom{T}, index::Index) where T
        new{T}(atom, AtomStatus.UNKNOWN, Set(), index)
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
end
function Line{T}() where T 
    return Line{T}(Vector{T}(undef, 3), Vector{T}(undef, 3))
end

struct Plane{T}
    r::Vector{T} # point 
    n::Vector{T} # normal
end
function Plane{T}(n::AbstractVector{T}, d::T) where T
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

function Plane{T}(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T}) where T
    return Plane{T}(a, cross(a.-b, b.-c))
end

struct ProbeSphere{T}
    r::Vector{T}
    radius::T
end

struct ProbePosition{T}
    status::Vector{ProbeStatus.T}
    points::Vector{Vector{T}}
end

mutable struct RSVertex{E<:AbstractEdge, F<:AbstractFace} <: AbstractVertex
    atom_index::Index
    edges::Set{E}
    faces::Set{F}

    index::Index
end

mutable struct RSEdge{T, V<:AbstractVertex, F<:AbstractFace} <: AbstractEdge
    vertices::NTuple{2, Union{Nothing, V}}
    faces::NTuple{2, Union{Nothing, F}}

    index::Index

    center_of_torus::Union{Nothing, Vector{T}}
    radius_of_torus::Union{Nothing, T}
    angle::Union{Nothing, T}
    circle0::Union{Nothing, Circle{T}}
    circle1::Union{Nothing, Circle{T}}
    intersection_point0::Union{Nothing, Vector{T}}
    intersection_point1::Union{Nothing, Vector{T}}
    singular::Union{Nothing, Bool}
end

mutable struct RSFace{T, V<:AbstractVertex, E<:AbstractEdge} <: AbstractFace
    vertices::NTuple{3, Union{Nothing, V}}
    edges::NTuple{3, Union{Nothing, E}}
    center::Vector{T}
    normal::Vector{T}
    singular::Bool

    index::Index
end

RSVertex(index) = RSVertex{RSEdge, RSFace}(index, Set(), Set(), 0)

RSVertex(v::RSVertex) = RSVertex{RSEdge, RSFace}(v.atom_index, v.edges, v.faces, v.index)

function RSEdge{T}() where T
    return RSEdge{T, RSVertex, RSFace}((nothing, nothing), (nothing, nothing), 0,
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end

function RSFace{T}() where T
    return RSFace{T, RSVertex, RSEdge}((nothing, nothing, nothing), (nothing, nothing, nothing), Vector{T}(), Vector{T}(), true, 0)
end

function RSFace(f::RSFace{T, RSVertex, RSEdge}) where T
    new{T, RSVertex, RSEdge}(f.vertices, f.edges, f.center, f.normal, f.singular, f.index)
end

function RSFace{T}(v1, v2, v3, e1, e2, e3, center, normal, singular, index) where T
    new{T, RSVertex, RSEdge}((v1, v2, v3), (e1, e2, e3), center, normal, singular, index)
end


# struct RSComputerStore
#     # Notes: Dicts can be handled as arrays with key==index. The key is the index of the atom in the large atom-array

#     # Properties that are related to a single atom can be stored directly in the atom
#     #atom_status::Dict{Index, AtomStatus.T}
#     #neighbours::Dict{Index, Set{Index}} # for each atom a list of its neighbours

#     #probe_radius::T
#     #neighboursOfTwoAtoms::Dict{Tuple{Index, Index}, Set{Index}} # for each pair of atoms (indices are ordered so that ind_a<ind_b) a list of their common neighbours
#     #probePositions::Dict{Tuple{Index, Index, Index}, Union{Nothing, ProbePosition}}

#     #newVertices::Set{RSVertex}
#     # newFaces::Set{RSFace}

#     vertices::Vector{Vector{RSVertex}}

#     #rm_vertices::Set{RSVertex}
# end

struct ReducedSurface{T}
    number_of_vertices::Index
    number_of_edges::Index
    number_of_faces::Index

    atoms::Vector{AtomWrapper{T}}

    vertices::Vector{RSVertex} #TODO this is not RSComputer->vertices_
    edges::Vector{Union{Nothing, RSEdge}}
    faces::Vector{Union{Nothing, RSFace}} 
end