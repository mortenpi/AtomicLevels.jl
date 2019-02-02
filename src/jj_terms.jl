function terms(orb::RelativisticOrbital, w::Int=one(Int))
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

# This is a workaround until seniority number are implemented for
# jj-coupled subshells.
intermediate_terms(orb::RelativisticOrbital, w::Int=one(Int)) =
    terms(orb, w)

export terms, intermediate_terms
