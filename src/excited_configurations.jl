function single_excitations!(excitations::Vector{Configuration{I,R}},
                             ref_set::Configuration{I,R},
                             orbitals::Vector{Orbital{I,R}},
                             min_occupancy::Vector{I},
                             max_occupancy::Vector{I},
                             excite_from::I,
                             keep_parity::Bool=true) where {I,R}
    for config ∈ excitations[end-excite_from+1:end]
        for (orb,occ,state) ∈ config
            state != :open && continue
            # If the orbital we propose to excite from is among those
            # from the reference set and already at its minimum
            # occupancy, we continue.
            i = findfirst(isequal(orb), ref_set.orbitals)
            !isnothing(i) && occ == min_occupancy[i] && continue
            for subs_orb ∈ orbitals
                subs_orb == orb && continue
                keep_parity && parity(subs_orb) != parity(orb) && continue
                # If the proposed substitution orbital is among those
                # from the reference set and already present in the
                # configuration to excite, we check if it is already
                # at its maximum occupancy, in which case we continue.
                j = findfirst(isequal(subs_orb), ref_set.orbitals)
                k = findfirst(isequal(subs_orb), config.orbitals)
                !isnothing(j) && !isnothing(k) && config.occupancy[k] == max_occupancy[j] && continue
                excited_config = replace(config, orb=>subs_orb)
                excited_config ∉ excitations && push!(excitations, excited_config)
            end
        end
    end
end

function excited_configurations(ref_set::Configuration{I,R},
                                orbitals::Orbital{I,R}...;
                                min_excitations::I=zero(I),
                                max_excitations::Union{I,Symbol}=:doubles,
                                min_occupancy::Vector{I}=zeros(I, length(peel(ref_set))),
                                max_occupancy::Vector{I}=[degeneracy(first(o)) for o in peel(ref_set)],
                                keep_parity::Bool=true) where {I<:Integer, R<:Rational{I}}

    if max_excitations isa Symbol
        max_excitations = if max_excitations == :singles
            1
        elseif max_excitations == :doubles
            2
        else
            throw(ArgumentError("Invalid maximum excitations specification $(max_excitations)"))
        end
    end

    all(min_occupancy .>= 0) || throw(ArgumentError("Invalid minimum occupancy requested"))
    all(max_occupancy .<= [degeneracy(first(o)) for o in peel(ref_set)]) ||
        throw(ArgumentError("Invalid maximum occupancy requested"))

    ref_set_core = core(ref_set)
    ref_set_peel = peel(ref_set)

    # The orbitals of the reference as valid orbitals to excite to as
    # well.
    orbitals = sort(unique(vcat(ref_set_peel.orbitals, orbitals...)))

    excitations = [ref_set_peel]
    excite_from = 1
    for i in 1:max_excitations
        single_excitations!(excitations, ref_set_peel, orbitals,
                            min_occupancy, max_occupancy, excite_from, keep_parity)
        excite_from = length(excitations)-excite_from
    end

    [ref_set_core + e for e in excitations]
end

export excited_configurations
