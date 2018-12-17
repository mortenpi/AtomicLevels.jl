# The main quantum number can either be an integer > 0 or a symbolic
# value (generally used to represent continuum electrons).
const MQ{I<:Integer} = Union{I,Symbol}

function check_orbital_ℓj(ℓ::I, j::R) where {I<:Integer, R<:Rational{I}}
    s = R(1,2)
    j₋ = abs(ℓ - s)
    j₊ = ℓ+s
    j ∈ j₋:j₊ || throw(ArgumentError("Invalid j = $(j) ∉ $(j₋):$(j₊)"))
end

struct Orbital{I<:Integer,N<:MQ{I}}
    n::N
    ℓ::I
    j::Rational{I}
    function Orbital(n::I, ℓ::I, j::R=ℓ+Rational{I}(1,2)) where {I<:Integer,R<:Rational{I}}
        n ≥ 1 || throw(ArgumentError("Invalid main quantum number $(n)"))
        0 ≤ ℓ && ℓ < n || throw(ArgumentError("Angular quantum number has to be ∈ [0,$(n-1)] when n = $(n)"))
        check_orbital_ℓj(ℓ, j)
        new{I,I}(n, ℓ, j)
    end
    function Orbital(n::Symbol, ℓ::I, j::R=ℓ+Rational{I}(1,2)) where {I<:Integer,R<:Rational{I}}
        check_orbital_ℓj(ℓ, j)
        new{I,Symbol}(n, ℓ, j)
    end
end

function Base.show(io::IO, orb::Orbital{I,N}) where {I,N}
    write(io, "$(orb.n)$(spectroscopic_label(orb.ℓ))")
    orb.j < orb.ℓ && write(io, "⁻")
end

flip_j(orb::Orbital) = Orbital(orb.n, orb.ℓ, orb.ℓ > 0 ? orb.ℓ - (orb.j - orb.ℓ) : orb.j)

degeneracy(orb::Orbital{I,N}) where {I,N} = 2orb.j + 1 |> I
non_rel_degeneracy(orb::Orbital) =
    orb.ℓ == 0 ? 2 : degeneracy(orb)+degeneracy(flip_j(orb))
parity(orb::Orbital{I,N}) where {I,N} = (-one(I))^orb.ℓ

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

parse_orbital_n(m::RegexMatch,i=1) =
    isnumeric(m[i][1]) ? parse(Int, m[i]) : Symbol(m[i])

function parse_orbital_ℓ(m::RegexMatch,i=2)
    ℓs = strip(m[i], ['[',']'])
    if isnumeric(ℓs[1])
        parse(Int, ℓs)
    else
        ℓi = findfirst(ℓs, spectroscopic)
        isnothing(ℓi) && throw(ArgumentError("Invalid spectroscopic label: $(m[i])"))
        first(ℓi) - 1
    end
end

function orbital_from_string(orb_str::AbstractString)
    m = match(r"^([0-9]+|.)([a-z]|\[[0-9]+\])([-]{0,1})$", orb_str)
    m === nothing && throw(ArgumentError("Invalid orbital string: $(orb_str)"))
    n = parse_orbital_n(m)
    ℓ = parse_orbital_ℓ(m)
    j = ℓ + (m[3] == "-" ? -1 : 1)*1//2
    Orbital(n, ℓ, j)
end

macro o_str(orb_str)
    orbital_from_string(orb_str)
end

function orbitals_from_string(orbs_str::AbstractString)
    map(split(orbs_str)) do orb_str
        m = match(r"^([0-9]+|.)\[([a-z]|[0-9]+)(-([a-z]|[0-9]+)){0,1}\]$", strip(orb_str))
        m === nothing && throw(ArgumentError("Invalid orbitals string: $(orb_str)"))
        n = parse_orbital_n(m)
        ℓs = map(filter(i -> !isnothing(m[i]), [2,4])) do i
            parse_orbital_ℓ(m, i)
        end
        orbs = map(ℓ -> Orbital(n, ℓ, ℓ-1//2), max(first(ℓs),1):last(ℓs))
        append!(orbs, map(ℓ -> Orbital(n, ℓ, ℓ+1//2), first(ℓs):last(ℓs)))
        sort(orbs)
    end |> o -> vcat(o...) |> sort
end

macro os_str(orbs_str)
    orbitals_from_string(orbs_str)
end

export Orbital, @o_str, @os_str, degeneracy, non_rel_degeneracy, parity
