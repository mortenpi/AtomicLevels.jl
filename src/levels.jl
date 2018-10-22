using UnicodeFun

struct Level
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
levels(conf::Config, term::Nothing) = [Level(conf,Term(0,0,1),0)]
levels(conf::Config) = sort(vcat([levels(conf,term) for term in terms(conf)]...))

import Base.show, Base.string

function latex(l::Level, rat_subscripts=true)
    J = if den(l.J) == 1
        "_{$(num(l.J))}"
    elseif rat_subscripts
        "_{$(num(l.J))/$(den(l.J))}"
    else
        "_{$(num(l.J))($(den(l.J)))}"
    end
    "|$(latex(l.conf))\\;$(latex(l.term))$(J)\\rangle"
end
show(io::IO, l::Level) = print(io, to_latex(latex(l, false)))

function show(io::IO, m::MIME"text/latex", l::Level, wrap = true)
    wrap && print(io, "\$")
    print(io, "\\mathrm{$(latex(l))}")
    wrap && print(io, "\$")
end

function string(l::Level)
    J = den(l.J) == 1 ? num(l.J) : replace("$(float(l.J))", ".", "_")
    "$(string(l.conf))_$(string(l.term))$(J)"
end

export Level, weight, ==, <, isless, hash, J_range, levels, show, string
