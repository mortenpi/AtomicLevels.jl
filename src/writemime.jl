function writemime{T<:Union{Term, Level}}(io::IO, m::MIME"text/latex", v::Vector{T})
    print(io, "\$[")
    for e in v[1:end-1]
        writemime(io, m, e, false)
        print(io, "\\quad")
    end
    writemime(io, m, v[end], false)
    print(io, "]\$")
end
