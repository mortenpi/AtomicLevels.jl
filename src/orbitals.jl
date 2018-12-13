# The main quantum number can either be an integer > 0 or a symbolic
# value (generally used to represent continuum electrons).
const MQ{I<:Integer} = Union{I,Symbol}

function check_orbital_ℓj(ℓ::I, j::R) where {I<:Integer, R<:Rational{I}}
    s = R(1,2)
    j₋ = abs(ℓ - s)
    j₊ = ℓ+s
    j ∈ j₋:j₊ || throw(ArgumentError("Invalid j = $(j) ∉ $(j₋):$(j₊)"))
end

struct Orbital{I<:Integer,R<:Rational{I},N<:MQ{I}}
    n::N
    ℓ::I
    j::R
    function Orbital(n::I, ℓ::I, j::R=ℓ+Rational{I}(1,2)) where {I<:Integer,R<:Rational{I}}
        n ≥ 1 || throw(ArgumentError("Invalid main quantum number $(n)"))
        0 ≤ ℓ && ℓ < n || throw(ArgumentError("Angular quantum number has to be ∈ [0,$(n-1)] when n = $(n)"))
        check_orbital_ℓj(ℓ, j)
        new{I,R,I}(n, ℓ, j)
    end
    function Orbital(n::Symbol, ℓ::I, j::R=ℓ+Rational{I}(1,2)) where {I<:Integer,R<:Rational{I}}
        check_orbital_ℓj(ℓ, j)
        new{I,R,Symbol}(n, ℓ, j)
    end
end

function Base.show(io::IO, orb::Orbital{I,R,N}) where {I,R,N}
    write(io, "$(orb.n)$(spectroscopic_label(orb.ℓ))")
    orb.j < orb.ℓ && write(io, "⁻")
end

flip_j(orb::Orbital) = Orbital(orb.n, orb.ℓ, orb.ℓ > 0 ? orb.ℓ - (orb.j - orb.ℓ) : orb.j)

degeneracy(orb::Orbital{I,R,N}) where {I,R,N} = 2orb.j + 1 |> I
non_rel_degeneracy(orb::Orbital) =
    orb.ℓ == 0 ? 2 : degeneracy(orb)+degeneracy(flip_j(orb))
parity(orb::Orbital{I,R,N}) where {I,R,N} = (-one(I))^orb.ℓ

nisless(an::T, bn::T) where T = an < bn
# Our convention is that symbolic quantum numbers are always greater
# than numeric ones, such that ks appears after 2p, etc.
nisless(an::I, bn::Symbol) where {I<:Integer} = true
nisless(an::Symbol, bn::I) where {I<:Integer} = false

function Base.isless(a::Orbital, b::Orbital)
    nisless(a.n, b.n) && return true
    a.n == b.n && a.ℓ < b.ℓ && return true
    a.n == b.n && a.ℓ == b.ℓ && a.j < b.j && return true
    false
end

function orbital_from_string(orb_str::AbstractString)
    m = match(r"^([0-9]+|.)([a-z]|\[[0-9]+\])([-]{0,1})$", orb_str)
    m === nothing && throw(ArgumentError("Invalid orbital string: $(orb_str)"))
    n = isnumeric(m[1][1]) ? parse(Int, m[1]) : Symbol(m[1])
    ℓi = findfirst(m[2], spectroscopic)
    ℓ = if ℓi === nothing
        parse(Int, strip(m[2], ['[', ']']))
    else
        first(ℓi) - 1
    end
    j = ℓ + (m[3] == "-" ? -1 : 1)*1//2
    Orbital(n, ℓ, j)
end

macro o_str(orb_str)
    orbital_from_string(orb_str)
end

export Orbital, @o_str, degeneracy, non_rel_degeneracy, parity
