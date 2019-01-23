using LinearAlgebra

@testset "jj2lsj" begin
    # These are not real tests, they only make sure that the
    # generation does not crash.
    R = jj2lsj(ros"1[s] 2[s-p] k[s-d]"...)
    # Since it is a rotation matrix, its inverse is simply the
    # transpose.
    @test norm(inv(Matrix(R)) - R') < 10eps()
end
