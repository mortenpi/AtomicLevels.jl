struct Configuration{O<:AbstractOrbital}
    orbitals::Vector{O}
    occupancy::Vector{Int}
    states::Vector{Symbol}
    function Configuration(
        orbitals::Vector{O},
        occupancy::Vector{Int},
        states::Vector{Symbol}=[:open for o in orbitals]) where {O<:AbstractOrbital}
        length(orbitals) == length(occupancy) ||
            throw(ArgumentError("Need to specify occupation numbers for all orbitals"))
        length(states) ≤ length(orbitals) ||
            throw(ArgumentError("Cannot provide more states than orbitals"))

        length(orbitals) == length(unique(orbitals)) ||
            throw(ArgumentError("Not all orbitals are unique"))

        if length(states) < length(orbitals)
            append!(states, repeat([:open], length(orbitals) - length(states)))
        end
        sd = setdiff(states, [:open, :closed, :inactive])
        isempty(sd) || throw(ArgumentError("Unknown orbital states $(sd)"))

        for i in eachindex(orbitals)
            occ = occupancy[i]
            occ < 1 && throw(ArgumentError("Invalid occupancy $(occ)"))
            orb = orbitals[i]
            degen_orb = degeneracy(orb)
            if occ > degen_orb
                if O <: Orbital || O <: SpinOrbital
                    throw(ArgumentError("Higher occupancy than possible for $(orb) with degeneracy $(degen_orb)"))
                elseif O <: RelativisticOrbital
                    orb_conj = flip_j(orb)
                    degen_orb_conj = degeneracy(orb_conj)
                    if occ == degen_orb + degen_orb_conj && orb_conj ∉ orbitals
                        occupancy[i] = degen_orb
                        push!(orbitals, orb_conj)
                        push!(occupancy, degen_orb_conj)
                        push!(states, states[i])
                    elseif orb_conj ∈ orbitals
                        throw(ArgumentError("Higher occupancy than possible for $(orb) with degeneracy $(degen_orb)"))
                    else
                        throw(ArgumentError("Can only specify higher occupancy for $(orb) if completely filling the $(orb),$(orb_conj) subshell (for total degeneracy of $(degen_orb+degen_orb_conj))"))
                    end
                end
            end
        end

        for i in eachindex(orbitals)
            states[i] == :closed && occupancy[i] < degeneracy(orbitals[i]) &&
                throw(ArgumentError("Can only close filled orbitals"))
        end

        p = sortperm(orbitals)

        new{O}(orbitals[p], occupancy[p], states[p])
    end
end

function Configuration(orbs::Vector{Tuple{O,Int,Symbol}}) where {O<:AbstractOrbital}
    orbitals = Vector{O}()
    occupancy = Vector{Int}()
    states = Vector{Symbol}()
    for (orb,occ,state) in orbs
        push!(orbitals, orb)
        push!(occupancy, occ)
        push!(states, state)
    end
    Configuration(orbitals, occupancy, states)
end

Configuration(orbital::O, occupancy::Int, state::Symbol=:open) where {O<:AbstractOrbital} =
    Configuration([orbital], [occupancy], [state])

Configuration{O}() where {O<:AbstractOrbital} =
    Configuration(O[], Int[], Symbol[])

const RelativisticConfiguration{N} = Configuration{RelativisticOrbital{N}}

issimilar(a::Configuration{<:O}, b::Configuration{<:O}) where {O<:AbstractOrbital} =
    a.orbitals == b.orbitals && a.occupancy == b.occupancy

Base.:(==)(a::Configuration{<:O}, b::Configuration{<:O}) where {O<:AbstractOrbital} =
    issimilar(a, b) && a.states == b.states

noble_gases = Dict(Orbital => Dict{String,Configuration{<:Orbital}}(),
                   RelativisticOrbital => Dict{String,Configuration{<:RelativisticOrbital}}())

"""
    fill(configuration)

Ensure all orbitals are at their maximum occupancy.
"""
function Base.fill(config::Configuration)
    map(config) do (orb,occ,state)
        orb,degeneracy(orb),state
    end |> Configuration
end

function Base.close(config::Configuration)
    map(config) do (orb,occ,state)
        orb,occ,:closed
    end |> Configuration
end

# This construct is needed since when showing configurations, they
# will be specialized on the Orbital parameterization, which we cannot
# index noble_gases with.
for O in [:Orbital,:RelativisticOrbital]
    @eval get_noble_gas(::Type{O}, k) where {O<:$O} = noble_gases[$O][k]
end

function write_orbitals(io::IO, config::Configuration)
    for (i,(orb,occ,state)) in enumerate(config)
        i > 1 && write(io, " ")
        write(io, "$(orb)")
        occ > 1 && write(io, to_superscript(occ))
        state == :closed && write(io, "ᶜ")
        state == :inactive && write(io, "ⁱ")
    end
end

function Base.show(io::IO, config::Configuration{O}) where O
    nc = length(config)
    if nc == 0
        write(io, "∅")
        return
    end
    core_config = core(config)
    ncc = length(core_config)
    if length(core_config) > 0
        core_printed = false
        for gas in ["Rn", "Xe", "Kr", "Ar", "Ne", "He"]
            gas_cfg = get_noble_gas(O, gas)
            ngc = length(gas_cfg)
            if ncc ≥ ngc && issimilar(core_config[1:length(gas_cfg)], gas_cfg)
                write(io, "[$(gas)]ᶜ")
                if ncc > ngc
                    write(io, " ")
                    write_orbitals(io, core_config[ngc+1:end])
                end
                nc > ncc && write(io, " ")
                core_printed = true
                break
            end
        end
        if !core_printed
            write_orbitals(io, core_config)
        end
    end
    write_orbitals(io, peel(config))
end

function state_sym(state::AbstractString)
    if state == "c"
        :closed
    elseif state == "i"
        :inactive
    else
        :open
    end
end

function core_configuration(::Type{O}, element::AbstractString, state::AbstractString) where {O<:AbstractOrbital}
    element ∉ keys(noble_gases[O]) && throw(ArgumentError("Unknown noble gas $(element)"))
    state = state_sym(state == "" ? "c" : state) # By default, we'd like cores to be frozen
    core_config = noble_gases[O][element]
    Configuration(core_config.orbitals, core_config.occupancy,
                  [state for o in core_config.orbitals])
end

function parse_orbital(::Type{O}, orb_str) where {O<:AbstractOrbital}
    m = match(r"^(([0-9]+|.)([a-z]|\[[0-9]+\])[-]{0,1})([0-9]*)([*ci]{0,1})$", orb_str)
    orbital_from_string(O, m[1]) , (m[4] == "") ? 1 : parse(Int, m[4]), state_sym(m[5])
end

function Base.parse(::Type{Configuration{O}}, conf_str::AbstractString) where {O<:AbstractOrbital}
    isempty(conf_str) && return Configuration{O}()
    orbs = split(conf_str, r"[\. ]")
    core_m = match(r"\[([a-zA-Z]+)\]([*ci]{0,1})", first(orbs))
    if core_m != nothing
        core_config = core_configuration(O, core_m[1], core_m[2])
        if length(orbs) > 1
            peel_config = Configuration(parse_orbital.(Ref(O), orbs[2:end]))
            Configuration(vcat(core_config.orbitals, peel_config.orbitals),
                          vcat(core_config.occupancy, peel_config.occupancy),
                          vcat(core_config.states, peel_config.states))
        else
            core_config
        end
    else
        Configuration(parse_orbital.(Ref(O), orbs))
    end
end

macro c_str(conf_str)
    parse(Configuration{Orbital}, conf_str)
end

macro rc_str(conf_str)
    parse(Configuration{RelativisticOrbital}, conf_str)
end

for O in [Orbital,RelativisticOrbital]
    for gas in ("He" => "1s2",
                "Ne" => "[He] 2s2 2p6",
                "Ar" => "[Ne] 3s2 3p6",
                "Kr" => "[Ar] 3d10 4s2 4p6",
                "Xe" => "[Kr] 4d10 5s2 5p6",
                "Rn" => "[Xe] 4f14 5d10 6s2 6p6")
        noble_gases[O][gas[1]] = parse(Configuration{O},gas[2])
    end
end

Base.getindex(conf::Configuration{O}, i::Integer) where O =
    (conf.orbitals[i], conf.occupancy[i], conf.states[i])
Base.getindex(conf::Configuration{O}, i::Union{<:UnitRange{<:Integer},<:AbstractVector{<:Integer}}) where O =
    Configuration([conf[ii] for ii in i])

Base.iterate(conf::Configuration{O}, (el, i)=(length(conf)>0 ? conf[1] : nothing,1)) where O =
    i > length(conf) ? nothing : (el, (conf[i==length(conf) ? i : i+1],i+1))

Base.length(conf::Configuration) = length(conf.orbitals)
Base.lastindex(conf::Configuration) = length(conf)
Base.eltype(conf::Configuration{O}) where O = (O,Int,Symbol)

function Base.isless(a::Configuration{<:O}, b::Configuration{<:O}) where {O<:AbstractOrbital}
    l = min(length(a),length(b))
    # If they are equal up to orbital l, designate the shorter config
    # as the smaller one.
    a[1:l] == b[1:l] && return length(a) == l
    norm_occ = (orb,w) -> 2w ≥ degeneracy(orb) ? degeneracy(orb) - w : w
    for ((orba,occa,statea),(orbb,occb,stateb)) in zip(a[1:l],b[1:l])
        if orba < orbb
            return true
        elseif orba == orbb
            # This is slightly arbitrary, but we designate the orbital
            # with occupation closest to a filled shell as the smaller
            # one.
            norm_occ(orba,occa) < norm_occ(orbb,occb) && return true
        else
            return false
        end
    end
    false
end

num_electrons(conf::Configuration) = sum(conf.occupancy)

Base.in(orb::Orbital, conf::Configuration{O}) where O =
    orb ∈ conf.orbitals

Base.filter(f::Function, conf::Configuration) =
    conf[filter(j -> f(conf[j]...), eachindex(conf.orbitals))]

core(conf::Configuration) = filter((orb,occ,state) -> state == :closed, conf)
peel(conf::Configuration) = filter((orb,occ,state) -> state != :closed, conf)
inactive(conf::Configuration) = filter((orb,occ,state) -> state == :inactive, conf)
active(conf::Configuration) = filter((orb,occ,state) -> state != :inactive, peel(conf))
bound(conf::Configuration) = filter((orb,occ,state) -> isbound(orb), conf)
continuum(conf::Configuration) = filter((orb,occ,state) -> !isbound(orb), peel(conf))

parity(conf::Configuration) = mapreduce(o -> parity(o[1])^o[2], *, conf)
Base.count(conf::Configuration) = mapreduce(o -> o[2], +, conf)

function Base.replace(conf::Configuration{O₁}, orbs::Pair{O₂,O₃}) where {O<:AbstractOrbital,O₁<:O,O₂<:O,O₃<:O}
    src,dest = orbs
    orbitals = promote_type(O₁,O₂,O₃)[]
    append!(orbitals, conf.orbitals)
    occupancy = copy(conf.occupancy)
    states = copy(conf.states)

    i = findfirst(isequal(src), orbitals)
    isnothing(i) && throw(ArgumentError("$(src) not present in $(conf)"))

    j = findfirst(isequal(dest), orbitals)
    if isnothing(j)
        push!(orbitals, dest)
        push!(occupancy, 1)
        push!(states, :open)
    else
        occupancy[j] == degeneracy(dest) &&
            throw(ArgumentError("$(dest) already maximally occupied in $(conf)"))
        occupancy[j] += 1
    end

    occupancy[i] -= 1
    if occupancy[i] == 0
        deleteat!(orbitals, i)
        deleteat!(occupancy, i)
        deleteat!(states, i)
    end

    Configuration(orbitals, occupancy, states)
end

function Base.:+(a::Configuration{O₁}, b::Configuration{O₂}) where {O<:AbstractOrbital,O₁<:O,O₂<:O}
    orbitals = promote_type(O₁,O₂)[]
    append!(orbitals, a.orbitals)
    occupancy = copy(a.occupancy)
    states = copy(a.states)
 for (orb,occ,state) in b
        i = findfirst(isequal(orb), orbitals)
        if isnothing(i)
            push!(orbitals, orb)
            push!(occupancy, occ)
            push!(states, state)
        else
            occupancy[i] += occ
            states[i] == state || throw(ArgumentError("Incompatible states for $(orb): $(states[i]) and $state"))
        end
    end
    Configuration(orbitals, occupancy, states)
end

"""
    ⊗(::Union{Configuration, Vector{Configuration}}, ::Union{Configuration, Vector{Configuration}})

Given two collections of `Configuration`s, it creates an array of `Configuration`s with all
possible juxtapositions of configurations from each collection.

# Examples

```jldoctest
julia> [c"1s", c"2s"] ⊗ [c"2p-", c"2p"]
4-element Array{Configuration{RelativisticOrbital},1}:
 1s 2p⁻
 1s 2p
 2s 2p⁻
 2s 2p

julia> c"1s" ⊗ [c"2s2", c"2s 2p-"]
2-element Array{Configuration{RelativisticOrbital},1}:
 1s 2s²
 1s 2s 2p⁻
```
"""
⊗(a::Vector{<:Configuration}, b::Vector{<:Configuration}) =
    [x+y for x in a for y in b]
⊗(a::Union{<:Configuration,Vector{<:Configuration}}, b::Configuration) =
    a ⊗ [b]
⊗(a::Configuration, b::Vector{<:Configuration}) =
    [a] ⊗ b

"""
    rconfigurations_from_orbital(n, ℓ, occupancy)

Generate all `Configuration`s with relativistic orbitals corresponding to the
non-relativistic orbital with `n` and `ℓ` quantum numbers, with given occupancy.

# Examples

```jldoctest
julia> rconfigurations_from_orbital(3, 1, 2)
3-element Array{Configuration{RelativisticOrbital},1}:
 3p⁻²
 3p⁻ 3p
 3p²
```
"""
function rconfigurations_from_orbital(n::N, ℓ::Int, occupancy::Int) where {N<:MQ}
    n isa Integer && ℓ + 1 > n && throw(ArgumentError("ℓ=$ℓ too high for given n=$n"))
    occupancy > 2*(2ℓ + 1) && throw(ArgumentError("occupancy=$occupancy too high for given ℓ=$ℓ"))

    degeneracy_ℓm, degeneracy_ℓp = 2ℓ, 2ℓ + 2 # degeneracies of nℓ- and nℓ orbitals
    nlow_min = max(occupancy - degeneracy_ℓp, 0)
    nlow_max = min(degeneracy_ℓm, occupancy)
    confs = RelativisticConfiguration[]
    for nlow = nlow_max:-1:nlow_min
        nhigh = occupancy - nlow
        conf = if nlow == 0
            Configuration([RelativisticOrbital(n, ℓ, ℓ + 1//2)], [nhigh])
        elseif nhigh == 0
            Configuration([RelativisticOrbital(n, ℓ, ℓ - 1//2)], [nlow])
        else
            Configuration(
                [RelativisticOrbital(n, ℓ, ℓ - 1//2),
                 RelativisticOrbital(n, ℓ, ℓ + 1//2)],
                [nlow, nhigh]
            )
        end
        push!(confs, conf)
    end
    return confs
end

"""
    rconfigurations_from_orbital(orbital::Orbital, occupancy)

Generate all `Configuration`s with relativistic orbitals corresponding to the
non-relativistic version of the `orbital` with a given occupancy.

# Examples

```jldoctest
julia> rconfigurations_from_orbital(o"3p", 2)
3-element Array{Configuration{RelativisticOrbital},1}:
 3p⁻²
 3p⁻ 3p
 3p²
```
"""
function rconfigurations_from_orbital(orbital::Orbital, occupation::Integer)
    rconfigurations_from_orbital(orbital.n, orbital.ℓ, occupation)
end

function rconfigurations_from_nrstring(orb_str::AbstractString)
    m = match(r"^([0-9]+|.)([a-z]+)([0-9]+)?$", orb_str)
    isnothing(m) && throw(ArgumentError("Invalid orbital string: $(orb_str)"))
    n = parse_orbital_n(m)
    ℓ = parse_orbital_ℓ(m)
    occupancy = isnothing(m[3]) ? 1 : parse(Int, m[3])
    return rconfigurations_from_orbital(n, ℓ, occupancy)
end

"""
    @rcs_str -> Vector{Configuration{RelativisticOrbital}}

Construct a `Vector` of all `Configuration`s corresponding to the non-relativistic `nℓ`
orbital with the given occupancy from the input string.

The string is assumed to have the following syntax: `\$(n)\$(ℓ)\$(occupancy)`, where `n`
and `occupancy` are integers, and `ℓ` is in spectroscopic notation.

# Examples

```jldoctest
julia> rcs"3p2"
3-element Array{Configuration{Int64},1}:
 3p⁻²
 3p⁻ 3p
 3p²
```
"""
macro rcs_str(s)
    rconfigurations_from_nrstring(s)
end


"""
    spin_configurations(configuration)

Generate all possible configurations of spin-orbitals from
`configuration`, i.e. all permissible values for the quantum numbers
`n`, `ℓ`, `mℓ`, `ms` for each electron. Example:

    spin_configuration(c"1s2") -> 1s₀α 1s₀β
"""
function spin_configurations(c::Configuration{<:Orbital})
    states = Dict{Orbital,Symbol}()
    orbitals = map(c) do (orb,occ,state)
        states[orb] = state
        sorbs = spin_orbitals(orb)
        collect(combinations(sorbs, occ)) |> Vector{Vector{SpinOrbital}}
    end
    map(allchoices(orbitals)) do choice
        c = vcat(choice...)
        s = [states[orb.orb] for orb in c]
        Configuration(c, ones(Int,length(c)), s)
    end
end

"""
    spin_configurations(configurations)

For each configuration in `configurations`, generate all possible
configurations of spin-orbitals.
"""
spin_configurations(cs::Vector{<:Configuration}) =
    sort(vcat(map(spin_configurations, cs)...))

function Base.show(io::IO, c::Configuration{<:SpinOrbital})
    orbitals = Dict{Orbital, Vector{<:SpinOrbital}}()
    core_orbitals = sort(unique(map(o -> o.orb, core(c).orbitals)))
    core_cfg = Configuration(core_orbitals, ones(Int, length(core_orbitals))) |> fill |> close

    if !isempty(core_cfg)
        show(io, core_cfg)
        write(io, " ")
    end
    for orb in peel(c).orbitals
        orbitals[orb.orb] = push!(get(orbitals, orb.orb, SpinOrbital[]), orb)
    end
    map(sort(collect(keys(orbitals)))) do orb
        ℓ = orb.ℓ
        g = degeneracy(orb)
        sub_shell = orbitals[orb]
        if length(sub_shell) == g
            format("{1:s}{2:s}", orb, to_superscript(g))
        else
            map(mℓrange(orb)) do mℓ
                mℓshell = findall(o -> o.mℓ == mℓ, sub_shell)
                if length(mℓshell) == 2
                    format("{1:s}{2:s}{3:s}", orb, to_subscript(mℓ), to_superscript(2))
                elseif length(mℓshell) == 1
                    string(sub_shell[mℓshell[1]])
                else
                    ""
                end
            end |> so -> join(filter(s -> !isempty(s), so), " ")
        end
    end |> so -> write(io, join(so, " "))
end

export Configuration, @c_str, @rc_str,
    num_electrons, core, peel, active, inactive, bound, continuum, parity, ⊗, @rcs_str,
    spin_configurations
