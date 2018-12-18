# Largely borrowed from WignerSymbols.jl (https://github.com/Jutho/WignerSymbols.jl)
# Copyright (c) 2017: Jutho Haegeman; MIT "Expat" License
"""
    struct HalfInteger <: Real

Represents half-integer values.
"""
struct HalfInteger <: Real
    num::Int
    function HalfInteger(; twoX = nothing)
        return new(twoX)
    end
end
HalfInteger(x::Real) = convert(HalfInteger, x)

Base.:+(a::HalfInteger, b::HalfInteger) = HalfInteger(twoX = a.num+b.num)
Base.:-(a::HalfInteger, b::HalfInteger) = HalfInteger(twoX = a.num-b.num)
Base.:-(a::HalfInteger) = HalfInteger(twoX = -a.num)
Base.:*(a::Integer, b::HalfInteger) = HalfInteger(twoX = a * b.num)
Base.:*(a::HalfInteger, b::Integer) = b * a
Base.:<=(a::HalfInteger, b::HalfInteger) = a.num <= b.num
Base.:<(a::HalfInteger, b::HalfInteger) = a.num < b.num
Base.one(::Type{HalfInteger}) = HalfInteger(twoX = 2)
Base.zero(::Type{HalfInteger}) = HalfInteger(twoX = 0)

Base.promote_rule(::Type{HalfInteger}, ::Type{<:Integer}) = HalfInteger
Base.promote_rule(::Type{HalfInteger}, T::Type{<:Rational}) = T
Base.promote_rule(::Type{HalfInteger}, T::Type{<:Real}) = T

Base.convert(::Type{HalfInteger}, n::Integer) = HalfInteger(twoX = 2*n)
function Base.convert(::Type{HalfInteger}, r::Rational)
    if r.den == 1
        return HalfInteger(twoX = 2*r.num)
    elseif r.den == 2
        return HalfInteger(twoX = r.num)
    else
        throw(InexactError(:HalfInteger, HalfInteger, r))
    end
end
function Base.convert(::Type{HalfInteger}, r::Real)
    num = 2*r
    if isinteger(num)
        return HalfInteger(twoX = convert(Int, num))
    else
        throw(InexactError(:HalfInteger, HalfInteger, r))
    end
end
Base.convert(T::Type{<:Integer}, s::HalfInteger) =
    iseven(s.num) ? convert(T, s.num>>1) : throw(InexactError(Symbol(T), T, s))
Base.convert(T::Type{<:Rational}, s::HalfInteger) = convert(T, s.num//2)
Base.convert(T::Type{<:Real}, s::HalfInteger) = convert(T, s.num/2)
Base.convert(::Type{HalfInteger}, s::HalfInteger) = s

Base.isinteger(a::HalfInteger) = iseven(a.num)
ishalfinteger(a::HalfInteger) = true
ishalfinteger(a::Integer) = true
ishalfinteger(a::Rational) = (a.den == 1) || (a.den == 2)
ishalfinteger(a::Real) = isinteger(2*a)

converthalfinteger(a::Number) = convert(HalfInteger, a)

Base.numerator(a::HalfInteger) = iseven(a.num) ? div(a.num, 2) : a.num
Base.denominator(a::HalfInteger) = iseven(a.num) ? 1 : 2

"""
    parse(HalfInteger, s)

Parses the string `s` into the corresponding `HalfInteger`-value. String can either be a
number or a fraction of the form `<x>/2`.
"""
function Base.parse(::Type{HalfInteger}, s::AbstractString)
    if in('/', s)
        J_numerator, J_denominator = split(s, '/'; limit=2)
        parse(Int, J_denominator) == 2 ||
            throw(ArgumentError("Denominator not 2 in HalfInteger string '$s'."))
        HalfInteger(parse(Int, J_numerator) // 2)
    elseif !isempty(strip(s))
        HalfInteger(parse(Int, s))
    else
        throw(ArgumentError("input string is empty or only contains whitespace"))
    end
end

macro hi_str(s)
    parse(HalfInteger, s)
end

Base.show(io::IO, am::HalfInteger) =
    print(io, iseven(am.num) ? "$(div(am.num, 2))" : "$(am.num)/2")

function Base.hash(a::HalfInteger, h::UInt)
    iseven(a.num) && return hash(a.num>>1, h)
    num, den = a.num, 2
    den = 1
    pow = -1
    if abs(num) < 9007199254740992
        return hash(ldexp(Float64(num),pow), h)
    end
    h = Base.hash_integer(den, h)
    h = Base.hash_integer(pow, h)
    h = Base.hash_integer(num, h)
    return h
end


struct HalfIntegerRange <: AbstractVector{HalfInteger}
    start :: HalfInteger
    stop :: HalfInteger

    function HalfIntegerRange(start::HalfInteger, stop::HalfInteger)
        (start <= stop) ||
            throw(ArgumentError("Second argument must be greater or equal to the first."))
        isinteger(stop - start) ||
            throw(ArgumentError("Two arguments must have integer difference."))
        return new(start, stop)
    end
end
Base.iterate(it::HalfIntegerRange) = (it.start, it.start + 1)
Base.iterate(it::HalfIntegerRange, s) = (s <= it.stop) ? (s, s+1) : nothing
Base.length(it::HalfIntegerRange) = convert(Int, it.stop - it.start) + 1
Base.size(it::HalfIntegerRange) = (length(it),)
function Base.getindex(it::HalfIntegerRange, i::Integer)
    1 <= i <= length(it) || throw(BoundsError(it, i))
    it.start + i - 1
end
Base.IteratorEltype(::HalfIntegerRange) = Base.HasEltype()
Base.eltype(::HalfIntegerRange) = HalfInteger
Base.IteratorSize(::HalfIntegerRange) = Base.HasLength()

Base.:(:)(i::HalfInteger, j::HalfInteger) = HalfIntegerRange(i, j)

export HalfInteger, @hi_str
