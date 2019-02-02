using AtomicLevels

function print_block(fun::Function,out=stdout)
    io = IOBuffer()
    fun(io)
    data = split(String(take!(io)), "\n")
    if length(data) == 1
        println(out, "[ $(data[1])")
    elseif length(data) > 1
        for (p,dl) in zip(vcat("⎡", repeat(["⎢"], length(data)-2), "⎣"),data)
            println(out, "$(p) $(dl)")
        end
    end
end

function print_states(cfgs::Vector{<:Configuration})
    foreach(cfgs) do cfg
        print_block() do io
            println(io, cfg)
            foreach(states.(csfs(cfg))) do ss
                print_block(io) do io
                    println(io, last(first(first(ss)).level.csf.terms))
                    foreach(ss) do s
                        print_block(io) do io
                            println(io)
                            foreach(sss -> println(io, sss), s)
                        end
                    end
                end
            end
        end
    end
end
print_states(cfgs::Configuration...) = print_states([cfgs...])

export print_states
