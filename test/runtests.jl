using MieScattering
using Test

custom_approx(a, b) = isapprox(a, b; rtol=1e-6, atol=1e-6)

@testset "MieScattering.jl" begin
    @testset "Bohren and Huffuman (1983)" begin
        res = bhmie(1.5, 1.5 + 0.01im)

        @test custom_approx(res.q_ext, 0.7949794469989591)
        @test custom_approx(res.q_sca, 0.7400000922059121)
        @test custom_approx(res.q_back, 0.12393910482931672)
        @test custom_approx(res.asymm, 0.5023608330153132)
    end
end
