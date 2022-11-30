using MieScattering
using Test

custom_approx(a, b) = isapprox(a, b; rtol=1e-5, atol=1e-5)

@testset "MieScattering.jl" begin
    @testset "Bohren and Huffuman (1983)" begin
        @testset "Homogeneous sphere" begin
            @testset "Appendix A" begin
                res = bhmie(1.5, 1.5 + 0.01im)

                @test custom_approx(res.q_ext, 0.7949794469989591)
                @test custom_approx(res.q_sca, 0.7400000922059121)
                @test custom_approx(res.q_back, 0.12393910482931672)
                @test custom_approx(res.asymm, 0.5023608330153132)
            end

            @testset "Custom" begin
                m = 1.55
                r = 0.525
                λ = 0.6328
                x = 2π * r / λ
                res = bhmie(x, m)

                @test custom_approx(res.q_ext, 3.10543)
                @test custom_approx(res.q_sca, 3.10543)
                @test custom_approx(res.q_back, 2.92534)
            end
        end

        @testset "Coated sphere" begin
            @testset "Appendix B" begin
                λ = 3.0
                r_core = 0.171
                r_mantle = 6.265
                x_core = 2π * r_core / λ
                x_mantle = 2π * r_mantle / λ
                m_core = 1.59 + 0.66im
                m_mantle = 1.409 + 0.1747im
                res = bhcoat(x_core, x_mantle, m_core, m_mantle)

                @test custom_approx(res.q_ext, 2.32803)
                @test custom_approx(res.q_sca, 1.14341)
                @test custom_approx(res.q_back, 0.0285099)
            end
        end
    end
end
