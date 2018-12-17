abstract type AbstractOrbital end

# The main quantum number can either be an integer > 0 or a symbolic
# value (generally used to represent continuum electrons).
const MQ = Union{Int,Symbol}

parity(orb::O) where {O<:AbstractOrbital} = p"odd"^orb.ℓ

nisless(an::T, bn::T) where T = an < bn
# Our convention is that symbolic main quantum numbers are always
# greater than numeric ones, such that ks appears after 2p, etc.
nisless(an::I, bn::Symbol) where {I<:Integer} = true
nisless(an::Symbol, bn::I) where {I<:Integer} = false

# * Non-relativistic orbital

struct Orbital{N<:MQ} <: AbstractOrbital
    n::N
    ℓ::Int
    function Orbital(n::Int, ℓ::Int)
        n ≥ 1 || throw(ArgumentError("Invalid main quantum number $(n)"))
        0 ≤ ℓ && ℓ < n || throw(ArgumentError("Angular quantum number has to be ∈ [0,$(n-1)] when n = $(n)"))
        new{Int}(n, ℓ)
    end
    function Orbital(n::Symbol, ℓ::Int)
        new{Symbol}(n, ℓ)
    end
end

Base.show(io::IO, orb::Orbital{N}) where N =
    write(io, "$(orb.n)$(spectroscopic_label(orb.ℓ))")

degeneracy(orb::Orbital) = 2*(2orb.ℓ + 1)

function Base.isless(a::Orbital, b::Orbital)
    nisless(a.n, b.n) && return true
    a.n == b.n && a.ℓ < b.ℓ && return true
    false
end

# * Relativistic orbital

function check_orbital_ℓj(ℓ::Int, j::R) where {R<:Rational{Int}}
    s = R(1,2)
    j₋ = abs(ℓ - s)
    j₊ = ℓ+s
    j ∈ j₋:j₊ || throw(ArgumentError("Invalid j = $(j) ∉ $(j₋):$(j₊)"))
end

struct RelativisticOrbital{N<:MQ} <: AbstractOrbital
    n::N
    ℓ::Int
    j::Rational{Int}
    function RelativisticOrbital(n::Int, ℓ::Int, j::R=ℓ+Rational{Int}(1,2)) where {R<:Rational{Int}}
        n ≥ 1 || throw(ArgumentError("Invalid main quantum number $(n)"))
        0 ≤ ℓ && ℓ < n || throw(ArgumentError("Angular quantum number has to be ∈ [0,$(n-1)] when n = $(n)"))
        check_orbital_ℓj(ℓ, j)
        new{Int}(n, ℓ, j)
    end
    function RelativisticOrbital(n::Symbol, ℓ::Int, j::R=ℓ+Rational{Int}(1,2)) where {R<:Rational{Int}}
        check_orbital_ℓj(ℓ, j)
        new{Symbol}(n, ℓ, j)
    end
end

function Base.show(io::IO, orb::RelativisticOrbital)
    write(io, "$(orb.n)$(spectroscopic_label(orb.ℓ))")
    orb.j < orb.ℓ && write(io, "⁻")
end

flip_j(orb::RelativisticOrbital) =
    RelativisticOrbital(orb.n, orb.ℓ, orb.ℓ > 0 ? orb.ℓ - (orb.j - orb.ℓ) : orb.j)

degeneracy(orb::RelativisticOrbital{N}) where N = 2orb.j + 1 |> Int

function Base.isless(a::RelativisticOrbital, b::RelativisticOrbital)
    nisless(a.n, b.n) && return true
    a.n == b.n && a.ℓ < b.ℓ && return true
    a.n == b.n && a.ℓ == b.ℓ && a.j < b.j && return true
    false
end

# * Orbital construction from strings

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

function orbital_from_string(::Type{O}, orb_str::AbstractString) where {O<:AbstractOrbital}
    m = match(r"^([0-9]+|.)([a-z]|\[[0-9]+\])([-]{0,1})$", orb_str)
    m === nothing && throw(ArgumentError("Invalid orbital string: $(orb_str)"))
    n = parse_orbital_n(m)
    ℓ = parse_orbital_ℓ(m)
    if O == RelativisticOrbital
        j = ℓ + (m[3] == "-" ? -1 : 1)*1//2
        O(n, ℓ, j)
    else
        m[3] == "" || throw(ArgumentError("Non-relativistic orbitals cannot have their spins explicitly specified"))
        O(n, ℓ)
    end
end

macro o_str(orb_str)
    orbital_from_string(Orbital, orb_str)
end

macro ro_str(orb_str)
    orbital_from_string(RelativisticOrbital, orb_str)
end

function orbitals_from_string(::Type{O}, orbs_str::AbstractString) where {O<:AbstractOrbital}
    map(split(orbs_str)) do orb_str
        m = match(r"^([0-9]+|.)\[([a-z]|[0-9]+)(-([a-z]|[0-9]+)){0,1}\]$", strip(orb_str))
        m === nothing && throw(ArgumentError("Invalid orbitals string: $(orb_str)"))
        n = parse_orbital_n(m)
        ℓs = map(filter(i -> !isnothing(m[i]), [2,4])) do i
            parse_orbital_ℓ(m, i)
        end
        orbs = if O == RelativisticOrbital
            orbs = map(ℓ -> O(n, ℓ, ℓ-1//2), max(first(ℓs),1):last(ℓs))
            append!(orbs, map(ℓ -> O(n, ℓ, ℓ+1//2), first(ℓs):last(ℓs)))
        else
            map(ℓ -> O(n, ℓ), first(ℓs):last(ℓs))
        end
        sort(orbs)
    end |> o -> vcat(o...) |> sort
end

macro os_str(orbs_str)
    orbitals_from_string(Orbital, orbs_str)
end

macro ros_str(orbs_str)
    orbitals_from_string(RelativisticOrbital, orbs_str)
end

export Orbital, RelativisticOrbital, @o_str, @ro_str, @os_str, @ros_str, degeneracy, parity
