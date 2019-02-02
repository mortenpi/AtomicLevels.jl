# Other utilities

## Parity

AtomicLevels defines the [`Parity`](@ref) type, which is used to represent the parity of
atomic states etc.

```@meta
DocTestSetup = quote
    using AtomicLevels
end
```

```@docs
Parity
@p_str
```

The parity values also define an algebra and an ordering:

```jldoctest
julia> p"odd" < p"even"
true

julia> p"even" * p"odd"
odd

julia> (p"odd")^3
odd

julia> -p"odd"
even
```

```@meta
DocTestSetup = nothing
```
