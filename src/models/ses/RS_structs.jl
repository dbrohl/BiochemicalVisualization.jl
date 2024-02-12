@enumx AtomStatus begin
    UNKNOWN
end

struct AtomWrapper{T}
    atom::Atom{T}
    status::AtomStatus
    p
end

struct Circle{T}
    r::Vector{T} # center
    n::Vector{T} # normal
    radius::T
end

struct Plane{T}
    p::Vector{T} # point #TODO StaticArrays?
    n::Vector{T} # normal

    function Plane{T}(n::Vector{T}, d::T) where T
        p = Vector{T}(undef, 3)
        if all(iszero_tolerance.(n))
            throw(DivideError())
        end
        if !iszero_tolerance(n[1])
            p = [-d / n[1], T(0), T(0)]
        elseif !iszero_tolerance(n[2])
            p = [T(0), -d / n[2], T(0)]
        elseif !iszero_tolerance(n[3])
            p = [T(0), T(0), -d / n[3]]
        end
        return Plane{T}(p, n)
    end
end