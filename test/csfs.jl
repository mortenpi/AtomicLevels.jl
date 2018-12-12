@testset "CSFs" begin
    @testset "Construction" begin
        # These are not real tests, they only make sure that the
        # generation does not crash.
        sort(vcat(csfs(c"3s 3p2")..., csfs(c"3s 3p- 3p")...))
        csfs(c"[Kr] 5s2 5p-2 5p3 6s")
    end
end
