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

!!! todo "TODO"

    Interface of `AbstractOrbital`?

The [`@c_str`](@ref) and [`@rc_str`](@ref) string macros can be used to conveniently
construct configurations:

```@docs
@c_str
@rc_str
```

## Interface

Various methods exist to query configurations for information.

```@docs
num_electrons(::Configuration)
```

```@meta
DocTestSetup = nothing
```
