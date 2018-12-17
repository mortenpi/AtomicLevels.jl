struct Parity
    p::Bool
end

macro p_str(ps)
    if ps == "even"
        Parity(true)
    elseif ps == "odd"
        Parity(false)
    else
        throw(ArgumentError("Invalid parity string $(ps)"))
    end
end

function Base.convert(::Type{Parity}, i::I) where {I<:Integer}
    i == 1 && return p"even"
    i == -1 && return p"odd"
    throw(ArgumentError("Don't know how to convert $(i) to parity"))
end

Base.iseven(p::Parity) = p.p
Base.isodd(p::Parity) = !p.p
Base.isless(a::Parity, b::Parity) = isodd(a) && iseven(b)

Base.:*(a::Parity, b::Parity) = Parity(a == b)
Base.:^(p::Parity, i::I) where {I<:Integer} =
    p.p || iseven(i) ? Parity(true) : Parity(false)
Base.:-(p::Parity) = Parity(!p.p)

Base.show(io::IO, p::Parity) =
    write(io, iseven(p) ? "even" : "odd")

UnicodeFun.to_superscript(p::Parity) = iseven(p) ? "" : "ᵒ"

export Parity, @p_str
