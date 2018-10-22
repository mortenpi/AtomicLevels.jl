import Base: show
function show(io::IO, m::MIME"text/latex", v::Vector{T}) where {T<:Union{Term, Level}}
    print(io, "\$[")
    for e in v[1:end-1]
        show(io, m, e, false)
        print(io, "\\quad")
    end
    show(io, m, v[end], false)
    print(io, "]\$")
end
