using SuperSub

type Level
    conf::Config
    term::Term
    J::Rational
end

weight(l::Level) = (2l.J + 1)

import Base.==
==(l1::Level, l2::Level) = ((l1.conf == l2.conf) && (l1.term == l2.term) && (l1.J == l2.J))

# It makes no sense to sort levels with different electron configurations
import Base.<
<(l1::Level, l2::Level) = ((l1.conf == l2.conf)) && (l1.term < l2.term) ||
    ((l1.conf == l2.conf)) && (l1.term == l2.term) && (l1.J < l2.J)
import Base.isless
isless(l1::Level, l2::Level) = (l1 < l2)

import Base.hash
hash(l::Level) = hash((hash(l.conf),hash(l.term),l.J))

J_range(term::Term) = abs(term.L-term.S):(term.L+term.S)

levels(conf::Config, term::Term) = sort([Level(conf,term,J) for J in J_range(term)])
levels(conf::Config, term::Void) = [Level(conf,Term(0,0,1),0)]
levels(conf::Config) = sort(vcat([levels(conf,term) for term in terms(conf)]...))

import Base.print, Base.show, Base.string

function print(io::IO, l::Level)
    print(io, "|", l.conf, " ", l.term)
    if den(l.J) == 1
        print(io, subscript(num(l.J)))
    else
        print(io, subscript(l.J))
    end
    print(io, "〉")
end
show(io::IO, l::Level) = print(io, l)

function show(io::IO, m::MIME"text/latex", l::Level, wrap = true)
    wrap && print(io, "\$")
    show(io, m, l.conf, false)
    print(io, "{\\;}")
    show(io, m, l.term, false)
    J = den(l.J) == 1 ? "$(num(l.J))" : "$(num(l.J))/$(den(l.J))"
    print(io, "_{$J}")
    wrap && print(io, "\$")
end

string(l::Level) = "$(string(l.conf))_$(string(l.term))$(l.J)"

export Level, weight, ==, <, isless, hash, J_range, levels, print, show, string
