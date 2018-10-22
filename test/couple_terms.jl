using AtomicLevels
using Test

import AtomicLevels: ells, ref_set_list

# Coupling two 2S terms should yields singlets and triplets of S,P,D
@test couple_terms(Term(1,1//2,1), Term(1,1//2,1)) ==
    sort([Term(2,0,1), Term(2,1,1), Term(0,0,1), Term(1,1,1), Term(1,0,1), Term(0,1,1)])


@test couple_terms(Term(0,1//2,1), Term(1,1//2,-1)) ==
    sort([Term(1,0,-1), Term(1,1,-1)])

function test_coupling(o1, o2)
    c1 = couple_terms(terms(ref_set_list(o1)), terms(ref_set_list(o2)))
    r = ref_set_list("$o1 $o2")
    c2 = terms(r)
    @test c1 == c2
end

test_coupling("1s", "2p")
test_coupling("2s", "2p")
test_coupling("2p", "2p")
test_coupling("2p", "3p")
test_coupling("2p", "3d")
