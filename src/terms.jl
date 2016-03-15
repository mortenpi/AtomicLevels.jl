using SuperSub

type Term
    L::Integer
    S::Rational
    parity::Integer
end

multiplicity(t::Term) = round(Int, 2t.S + 1)
weight(t::Term) = (2t.L + 1) * multiplicity(t)

import Base.==
==(t1::Term, t2::Term) = ((t1.L == t2.L) && (t1.S == t2.S) && (t1.parity == t2.parity))

import Base.<
<(t1::Term, t2::Term) = ((t1.L < t2.L) || (t1.S < t2.S) || (t1.parity < t2.parity))
import Base.isless
isless(t1::Term, t2::Term) = (t1 < t2)

import Base.hash
hash(t::Term) = hash((t.L,t.S,t.parity))

function couple_terms(t1::Term, t2::Term)
    L1 = t1.L
    L2 = t2.L
    S1 = t1.S
    S2 = t2.S
    p = t1.parity * t2.parity
    vcat([[Term(L, S, p) for S in abs(S1-S2):(S1+S2)]
          for L in abs(L1-L2):(L1+L2)]...)
end

function couple_terms(t1s::Vector{Term}, t2s::Vector{Term})
    ts = map(t1s) do t1
        map(t2s) do t2
            couple_terms(t1, t2)
        end
    end
    sort(unique(vcat(vcat(ts...)...)))
end

couple_terms(ts::Vector{Vector{Term}}) = reduce(couple_terms, ts)

function couple_terms(ts::Vector{Void})
    [nothing]
end

function terms(orb::Orbital, p::Integer)
    ell = orb[2]
    occ = orb[3]
    g = degeneracy(orb)
    (occ > g/2) && (occ = g - occ)

    if occ == 1
        return [Term(ell,1//2,p)]
    elseif occ == 2
        SL = filter!(sl -> sum(sl)%2 == 0,
                     vcat([[(S,L) for S in [0, 1]] for L in 0:2ell]...))
        return sort([Term(L, S, p) for (S,L) in SL])
    end

    warn("No implementation for more than 2 equivalent electrons")
end

function terms(config::Config)
    config = filter(o -> o[3]<degeneracy(o), config)

    # All subshells are filled, total angular momentum = 0, total spin = 0
    if length(config) == 0
        return [Term(0, 0, 1)]
    end

    p = parity(config)

    ts = map(config) do orb
        terms(orb, p)
    end

    couple_terms(ts)
end
terms(config::AbstractString) = terms(ref_set_list(config))

import Base.print, Base.show, Base.writemime, Base.string

function ELL(L::Integer)
    if L<length(ells)
        uppercase(ells[L+1])
    else
        string(L)
    end
end

function print(io::IO, t::Term)
    print(io, superscript(multiplicity(t)))
    print(io, ELL(t.L))
    t.parity == -1  && print(io, "⁻")
end
show(io::IO, t::Term) = print(io, t)

function writemime(io::IO, ::MIME"text/latex", t::Term, wrap = true)
    wrap && print(io, "\$")
    p = t.parity == -1 ? "^{\\mathrm{o}}" : ""
    print(io, "^{$(multiplicity(t))}\\mathrm{$(ELL(t.L))}$p")
    wrap && print(io, "\$")
end

string(t::Term) = "$(2t.S+1)$(ELL(t.L))$(t.parity == -1 ? '-' : "")"

export Term, multiplicity, weight, ==, <, isless, hash, couple_terms, terms, print, show, string
