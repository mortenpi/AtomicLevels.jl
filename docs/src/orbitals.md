# Atomic orbitals

```@meta
DocTestSetup = quote
    using AtomicLevels
end
```

Atomic orbitals, i.e. single particle states with well-defined orbital or total angular
momentum, are usually the basic building blocks of atomic states. AtomicLevels provide
various types and methods to work with orbitals.

## Orbital types

AtomicLevels provides two basic types for labelling atomic orbitals: [`Orbital`](@ref) and
[`RelativisticOrbital`](@ref).

```@docs
Orbital
RelativisticOrbital
```

To have labels for orbitals with all quantum numbers specified (i.e. including ``m_\ell``
and ``m_s``), the [`SpinOrbital`](@ref) type can be used.

```@docs
SpinOrbital
```

The string macros [`@o_str`](@ref) and [`@ro_str`](@ref) can be used to conveniently contruct
orbitals, while [`@os_str`](@ref) and [`@ros_str`](@ref) can be used to construct whole lists
of them very easily.

```@docs
@o_str
@ro_str
@os_str
@ros_str
```

## Methods

```@docs
isless
degeneracy
parity(::Orbital)
symmetry
isbound
mℓrange
spin_orbitals
@κ_str
```

```@meta
DocTestSetup = nothing
```
