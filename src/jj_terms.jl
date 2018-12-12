using Parameters

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

function jj_terms(orb::Orbital{I,R}, w::I=one(I)) where {I,R}
    @unpack n,ℓ,j = orb

    2w ≥ 2j+1 && (w = 2j+1-w)
    w == 0 && return [zero(R)]
    w == 1 && return [j]

    # Forms full Cartisian product of all mⱼ, not necesserily the most
    # performant method.
    mⱼs = filter(allunique, collect(allchoices([-j:j for i = 1:w])))
    MJs = map(sum, mⱼs)

    Js = R[]

    while !isempty(MJs)
        # Identify the maximum MJ and associate it with J.
        MJmax = maximum(MJs)
        n = count(isequal(MJmax), MJs)
        append!(Js, repeat([MJmax], n))
        # Remove the -MJ:MJ series, n times.
        for MJ = -MJmax:MJmax
            deleteat!(MJs, findall(isequal(MJ), MJs)[1:n])
        end
    end

    # Do we really want unique here?
    sort(unique(Js))
end

export jj_terms
