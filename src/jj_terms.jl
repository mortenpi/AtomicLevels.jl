couple_terms(J1::T, J2::T) where {T <: Union{Integer,HalfInteger}} =
    collect(abs(J1-J2):(J1+J2))
function couple_terms(J::Vector{T}, j₀::T=zero(T)) where {T <: Union{Integer,HalfInteger}}
    ts = Vector{Vector{T}}()
    for t in couple_terms(j₀, J[1])
        if length(J) == 1
            push!(ts, [j₀, t])
        else
            for ts′ in couple_terms(J[2:end], t)
                push!(ts, vcat(j₀, ts′...))
            end
        end
    end
    ts
end

couple_terms(J1::Real, J2::Real) =
    couple_terms(convert(HalfInteger, J1), convert(HalfInteger, J2))
couple_terms(J::Vector{T}, j₀::Real=zero(T)) where {T <: Real} =
    couple_terms(convert.(HalfInteger, J), convert(HalfInteger, j₀))

function jj_terms(orb::RelativisticOrbital, w::Int=one(Int))
    j = kappa_to_j(orb.κ)
    w <= 2j+1 || throw(ArgumentError("w=$w too large for $orb orbital"))

    2w ≥ 2j+1 && (w = convert(Int, 2j) + 1 - w)
    w == 0 && return [zero(HalfInteger)]
    w == 1 && return [j]

    # Forms full Cartesian product of all mⱼ, not necessarily the most
    # performant method.
    mⱼs = filter(allunique, collect(allchoices([-j:j for i = 1:w])))
    MJs = map(x -> reduce(+, x), mⱼs) # TODO: make map(sum, mⱼs) work

    Js = HalfInteger[]

    while !isempty(MJs)
        # Identify the maximum MJ and associate it with J.
        MJmax = maximum(MJs)
        N = count(isequal(MJmax), MJs)
        append!(Js, repeat([MJmax], N))
        # Remove the -MJ:MJ series, N times.
        for MJ = -MJmax:MJmax
            deleteat!(MJs, findall(isequal(MJ), MJs)[1:N])
        end
    end

    # Do we really want unique here?
    sort(unique(Js))
end

export jj_terms
