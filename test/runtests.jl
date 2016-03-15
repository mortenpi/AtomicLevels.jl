using AtomicLevels
using Base.Test

include("gen_terms.jl")
# Coupling two 2S terms should yields singlets and triplets of S,P,D
@test couple_terms(Term(1,1//2,1), Term(1,1//2,1)) == sort([Term(2,0,1), Term(2,1,1), Term(0,0,1), Term(1,1,1), Term(1,0,1), Term(0,1,1)])
