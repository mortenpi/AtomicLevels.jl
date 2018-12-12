struct allchoices{T,V<:AbstractVector{T}}
    choices::Vector{V}
    dims::Vector{Int}
end
function allchoices(choices::Vector{V}) where {T,V<:AbstractVector{T}}
    allchoices(choices, length.(choices))
end
Base.length(ac::allchoices) = prod(ac.dims)

function Base.iterate(ac::allchoices{T,V}, (el,i)=(first.(ac.choices),ones(Int,length(ac.dims)))) where {T,V}
    i == 0 && return nothing
    finish_next = false
    for j ∈ reverse(eachindex(i))
        i[j] += 1
        if i[j] > ac.dims[j]
            j == 1 && (finish_next = true)
            i[j] = 1
        else
            break
        end
    end
    el, ([c[ii] for (c,ii) in zip(ac.choices,i)], finish_next ? 0 : i)
end
