using SuperSub

type Term
    L::Integer
    S::Rational
    parity::Integer
end

function term_string(s::AbstractString)
    T_pat = r"([0-9]+)([A-Z])([oe ]{0,1})"
    ismatch(T_pat, s) || error("Invalid term string, $s")
    m = match(T_pat, s)
    L = lowercase(m[2][1])
    Term(findfirst(ells, L)-1, (parse(Int, m[1]) - 1)//2,
         m[3] == "o" ? -1 : 1)
end

macro T_str(s::AbstractString)
    term_string(s)
end

multiplicity(t::Term) = round(Int, 2t.S + 1)
weight(t::Term) = (2t.L + 1) * multiplicity(t)

import Base.==
==(t1::Term, t2::Term) = ((t1.L == t2.L) && (t1.S == t2.S) && (t1.parity == t2.parity))

import Base.<
<(t1::Term, t2::Term) = ((t1.S < t2.S) || (t1.S == t2.S) && (t1.L < t2.L)
                         || (t1.S == t2.S) && (t1.L == t2.L) && (t1.parity < t2.parity))
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
    sort(vcat([[Term(L, S, p) for S in abs(S1-S2):(S1+S2)]
               for L in abs(L1-L2):(L1+L2)]...))
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

include("xu2006.jl")

# This function calculates the term symbol for a given orbital ℓʷ
function xu_terms(ell::Integer, w::Integer, p::Integer)
    ts = map(((w//2 - floor(Int, w//2)):w//2)) do S
        S_p = round(Int,2S)
        map(0:w*ell) do L
            repmat([Term(L,S,p)], Xu.X(w,ell, S_p, L))
        end
    end
    vcat(vcat(ts...)...)
end

function terms(orb::Orbital)
    ell = orb[2]
    w = orb[3]
    g = degeneracy(orb)
    (w > g/2) && (w = g - w)

    p = parity(orb)

    if w == 1
        return [Term(ell,1//2,p)]
    end

    xu_terms(ell, w, p)
end

function terms(config::Config)
    config = filter(o -> o[3]<degeneracy(o), config)

    # All subshells are filled, total angular momentum = 0, total spin = 0
    if length(config) == 0
        return [Term(0, 0, 1)]
    end

    ts = map(config) do orb
        terms(orb)
    end

    couple_terms(ts)
end
terms(config::AbstractString) = terms(ref_set_list(config))

import Base.print, Base.show, Base.string

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
    t.parity == -1  && print(io, "ᵒ")
end
show(io::IO, t::Term) = print(io, t)

function show(io::IO, ::MIME"text/latex", t::Term, wrap = true)
    wrap && print(io, "\$")
    p = t.parity == -1 ? "^{\\mathrm{o}}" : ""
    print(io, "^{$(multiplicity(t))}\\mathrm{$(ELL(t.L))}$p")
    wrap && print(io, "\$")
end

string(t::Term, show_parity::Bool = true) = "$(round(Int, 2t.S+1))$(ELL(t.L))$(show_parity ? (t.parity == -1 ? 'o' : 'e') : "")"

export Term, term_string, @T_str, multiplicity, weight, ==, <, isless, hash, couple_terms, terms, print, show, string
