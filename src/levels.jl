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
<(l1::Level, l2::Level) = ((l1.conf == l2.conf)) && (l1.term < l2.term) || (l1.J < l2.J)
import Base.isless
isless(l1::Level, l2::Level) = (l1 < l2)

import Base.hash
hash(l::Level) = hash((hash(l.config),hash(l.term),l.J))

function levels(conf::Config, term::Term)
    L = term.L
    S = term.S
    sort([Level(conf,term,J) for J in abs(L-S):(L+S)])
end

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

string(l::Level) = "$(string(l.conf))_$(string(l.term))$(l.J)"

export Level, weight, ==, <, isless, hash, levels, print, show, string
