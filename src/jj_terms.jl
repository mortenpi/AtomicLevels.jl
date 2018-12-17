couple_terms(J1::N, J2::N) where {I<:Integer,R<:Rational{I},N<:Union{I,R}} =
    collect(abs(J1-J2):(J1+J2))

function couple_terms(J::Vector{N}, j₀::N=zero(N)) where {I<:Integer,R<:Rational{I},N<:Union{I,R}}
    ts = Vector{Vector{N}}()
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

function jj_terms(orb::RelativisticOrbital, w::Int=one(Int))
    @unpack ℓ,j = orb

    2w ≥ 2j+1 && (w = 2j+1-w)
    w == 0 && return [zero(Rational{Int})]
    w == 1 && return [j]

    # Forms full Cartesian product of all mⱼ, not necessarily the most
    # performant method.
    mⱼs = filter(allunique, collect(allchoices([-j:j for i = 1:w])))
    MJs = map(sum, mⱼs)

    Js = Rational{Int}[]

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
