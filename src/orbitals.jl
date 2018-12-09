struct Orbital{I<:Integer,R<:Rational{I}}
    n::I
    ℓ::I
    j::R
    function Orbital(n::I, ℓ::I, j::R=ℓ+Rational{I}(1,2)) where {I<:Integer,R<:Rational{I}}
        n ≥ 1 || throw(ArgumentError("Invalid main quantum number $(n)"))
        0 ≤ ℓ && ℓ < n || throw(ArgumentError("Angular quantum number has to be ∈ [0,$(n-1)] when n = $(n)"))
        s = R(1,2)
        j₋ = abs(ℓ - s)
        j₊ = ℓ+s
        j ∈ j₋:j₊ || throw(ArgumentError("Invalid j = $(j) ∉ $(j₋):$(j₊)"))
        new{I,R}(n, ℓ, j)
    end
end

function Base.show(io::IO, orb::Orbital{I,R}) where {I,R}
    write(io, "$(orb.n)$(spectroscopic_label(orb.ℓ))")
    orb.j < orb.ℓ && write(io, "⁻")
end

degeneracy(orb::Orbital{I,R}) where {I,R} = 2orb.j + 1 |> I
parity(orb::Orbital{I,R}) where {I,R} = (-one(I))^orb.ℓ

flip_j(orb::Orbital) = Orbital(orb.n, orb.ℓ, orb.ℓ > 0 ? orb.ℓ - (orb.j - orb.ℓ) : orb.j)

function Base.isless(a::Orbital, b::Orbital)
    a.n < b.n && return true
    a.n == b.n  && a.ℓ < b.ℓ && return true
    a.n == b.n  && a.ℓ == b.ℓ && a.j < b.j && return true
    false
end

function orbital_from_string(orb_str::AbstractString)
    m = match(r"([0-9]+)([a-z]|\[[0-9]+\])([-]{0,1})", orb_str)
    m === nothing && throw(ArgumentError("Invalid orbital string: $(orb_str)"))
    n = parse(Int, m[1])
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

export Orbital, @o_str, degeneracy, parity
