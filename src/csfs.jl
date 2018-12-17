struct CSF{I<:Integer, T<:Union{Term{I,Rational{I}},Rational{I}}}
    config::Configuration{I}
    subshell_terms::Vector{T}
    terms::Vector{T}

    function CSF(config::Configuration, subshell_terms::Vector{<:Term}, terms::Vector{<:Term})
        error("Non-relativistic CSFs not yet implemented")
    end

    function CSF(config::Configuration{I}, subshell_terms::Vector{R}, terms::Vector{R}) where {I,R<:Rational{I}}
        length(subshell_terms) == length(peel(config)) ||
            throw(ArgumentError("Need to provide $(length(peel(config))) subshell terms for $(config)"))
        length(terms) == length(peel(config)) ||
            throw(ArgumentError("Need to provide $(length(peel(config))) terms for $(config)"))
        new{I,R}(config, subshell_terms, terms)
    end
end

# We possibly want to sort on configuration as well
Base.isless(a::CSF, b::CSF) = last(a.terms) < last(b.terms)

function csfs(config::Configuration{I}) where I
    p_config = peel(config)
    js = map(p_config) do (orb,occ,state)
        jj_terms(orb,occ)
    end
    map(allchoices(js)) do j
        map(couple_terms(j)) do term
            CSF(config, j, term[2:end])
        end
    end |> c -> vcat(c...) |> sort
end

function Base.show(io::IO, csf::CSF{I,R}) where {I<:Integer,R<:Rational{I}}
    c = core(csf.config)
    cw = length(string(c))

    w = 0
    p = peel(csf.config)
    for (orb,occ,state) in p
        w = max(w, length(string(orb))+length(digits(occ)))
    end

    w += 1

    cfmt = FormatExpr("{1:$(cw)s} ")
    ofmt = FormatExpr("{1:<$(w+1)s}")
    ifmt = FormatExpr("{1:>$(w+1)d}")
    rfmt = FormatExpr("{1:>$(w-1)d}/2")

    printfmt(io, cfmt, c)
    for (orb,occ,state) in p
        printfmt(io, ofmt, join([string(orb), occ > 1 ? to_superscript(occ) : ""], ""))
    end
    println(io)

    printfmt(io, cfmt, "")
    for j in csf.subshell_terms
        if denominator(j) == 1
            printfmt(io, ifmt, numerator(j))
        else
            printfmt(io, rfmt, numerator(j))
        end
    end
    println(io)

    printfmt(io, cfmt, "")
    print(io, " ")
    for j in csf.terms
        if denominator(j) == 1
            printfmt(io, ifmt, numerator(j))
        else
            printfmt(io, rfmt, numerator(j))
        end
    end
    print(io, parity(csf.config) == 1 ? "+" : "-")
end

export CSF, csfs
