struct Term
    L::HalfInteger
    S::HalfInteger
    parity::Parity
end
Term(L::Real, S::Real, parity::Integer) =
    Term(convert(HalfInteger, L), convert(HalfInteger, S), convert(Parity, parity))

function Base.parse(::Type{Term}, s::AbstractString)
    m = match(r"([0-9]+)([A-Z]|\[[0-9/]+\])([oe ]{0,1})", s)
    isnothing(m) && throw(ArgumentError("Invalid term string $s"))
    L = lowercase(m[2])
    L = if L[1] == '['
        L = strip(L, ['[',']'])
        if occursin("/", L)
            Ls = split(L, "/")
            length(Ls) == 2 && length(Ls[1]) > 0 && length(Ls[2]) > 0 ||
                throw(ArgumentError("Invalid term string $(s)"))
            Rational(parse(Int, Ls[1]),parse(Int, Ls[2]))
        else
            parse(Int, L)
        end
    else
        findfirst(L, spectroscopic)[1]-1
    end
    denominator(L) ∈ [1,2] || throw(ArgumentError("L must be integer or half-integer"))
    S = (parse(Int, m[1]) - 1)//2
    Term(L, S, m[3] == "o" ? p"odd" : p"even")
end

macro T_str(s::AbstractString)
    parse(Term, s)
end

multiplicity(t::Term) = convert(Int, 2t.S + 1)
weight(t::Term) = (2t.L + 1) * multiplicity(t)

import Base.==
==(t1::Term, t2::Term) = ((t1.L == t2.L) && (t1.S == t2.S) && (t1.parity == t2.parity))

import Base.<
<(t1::Term, t2::Term) = ((t1.S < t2.S) || (t1.S == t2.S) && (t1.L < t2.L)
                         || (t1.S == t2.S) && (t1.L == t2.L) && (t1.parity < t2.parity))
import Base.isless
isless(t1::Term, t2::Term) = (t1 < t2)

Base.hash(t::Term) = hash((t.L,t.S,t.parity))

function couple_terms(t1::Term, t2::Term)
    # It is assumed that t1 and t2 originate from non-equivalent
    # electrons, since the vector model does not predict correct term
    # couplings for equivalent electrons; some of the generated terms
    # would violate the Pauli principle; cf. Cowan p. 108–109.
    L1 = t1.L
    L2 = t2.L
    S1 = t1.S
    S2 = t2.S
    p = t1.parity * t2.parity
    sort(vcat([[Term(L, S, p) for S in abs(S1-S2):(S1+S2)]
               for L in abs(L1-L2):(L1+L2)]...))
end

function couple_terms(t1s::Vector{<:Term}, t2s::Vector{<:Term})
    ts = map(t1s) do t1
        map(t2s) do t2
            couple_terms(t1, t2)
        end
    end
    sort(unique(vcat(vcat(ts...)...)))
end

couple_terms(ts::Vector{<:Vector{<:Term}}) = foldl(couple_terms, ts)

include("xu2006.jl")

# This function calculates the term symbol for a given orbital ℓʷ
function xu_terms(ℓ::Int, w::Int, p::Parity)
    ts = map(((w//2 - floor(Int, w//2)):w//2)) do S
        S′ = 2S |> Int
        map(L -> repeat([Term(L,S,p)], Xu.X(w,ℓ, S′, L)), 0:w*ℓ)
    end
    vcat(vcat(ts...)...)
end

function terms(orb::Orbital, occ::Int)
    ℓ = orb.ℓ
    g = degeneracy(orb)
    occ > g && throw(ArgumentError("Invalid occupancy $occ for $orb with degeneracy $g"))
    (occ > g/2 && occ != g) && (occ = g - occ)

    p = parity(orb)^occ
    if occ == 1
        [Term(ℓ,1//2,p)] # Single electron
    elseif ℓ == 0 && occ == 2 || occ == degeneracy(orb)
        [Term(0,0,p"even")] # Filled ℓ shell
    else
        xu_terms(ℓ, occ, p) # All other cases
    end
end

function terms(config::Configuration{O}) where {O<:Orbital}
    ts = map(config) do (orb,occ,state)
        terms(orb,occ)
    end

    couple_terms(ts)
end

function write_L(io::IO, term::Term)
    if isinteger(term.L)
        write(io, uppercase(spectroscopic_label(convert(Int, term.L))))
    else
        write(io, "[$(numerator(term.L))/$(denominator(term.L))]")
    end
end

function Base.show(io::IO, term::Term)
    write(io, to_superscript(multiplicity(term)))
    write_L(io, term)
    write(io, to_superscript(term.parity))
end

export Term, @T_str, multiplicity, weight, couple_terms, terms
