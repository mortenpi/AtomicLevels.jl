# Atomic configurations

```@meta
DocTestSetup = quote
    using AtomicLevels
end
```

We define a configuration to be a set of orbitals with their associated occupation (i.e.
the number of electron on that orbital).
We can represent a particular configuration with an instance of the [`Configuration`](@ref)
type.

```@docs
Configuration
```

The [`@c_str`](@ref) and [`@rc_str`](@ref) string macros can be used to conveniently
construct configurations:

```@docs
@c_str
@rc_str
```

## Interface

Various methods exist to manipulate or transform configurations or to query them for
information.

```@docs
num_electrons(::Configuration)
Base.delete!
Base.:(-)
Base.close
close!
Base.fill
Base.fill!
```

## Generating configuration lists

The [`⊗`](@ref) operator can be used to easily generate lists of configurations from existing
pieces. E.g. to create all the valence configurations on top of an closed core, you only
need to write

```jldoctest
julia> c"[Ne]" ⊗ [c"3s2", c"3s 3p", c"3p2"]
3-element Array{Configuration{Orbital{Int64}},1}:
 [Ne]ᶜ 3s²
 [Ne]ᶜ 3s 3p
 [Ne]ᶜ 3p²
```

That can be combined with the [`@rcs_str`](@ref) string macro to easily generate all possible
relativistic configurations from a non-relativistic definition:

```jldoctest
julia> rc"[Ne] 3s2" ⊗ rcs"3p2"
3-element Array{Configuration{RelativisticOrbital{Int64}},1}:
 [Ne]ᶜ 3s² 3p⁻²
 [Ne]ᶜ 3s² 3p⁻ 3p
 [Ne]ᶜ 3s² 3p²
```

```@docs
⊗
@rcs_str
```

## Spin configurations

```@docs
spin_configurations
substitutions
```

```@meta
DocTestSetup = nothing
```
