export bhmie

@doc raw"""
`bhmie([T=Float64,], x::Real, m::Number, Nθ::Integer=181; nextra::Integer=15, custom_nstop=0)`

Inputs:

- `T`: Type used for calculation. All real numbers will be stored as `T`, while all complex numbers will be stored as `Complex{T}`.
- `x`: Size parameter of the sphere scatterer. Defined as ``\frac{2\pi r}{\lambda}``
- `m`: Relative refractive index of the scatterer.
- `Nθ`: Number of scattering angles to calculate. Default is `181`.

Additional keyword arguments:

- `nextra`: Extra terms used for the downward calculation of the `d` function. Default is `15`.
- `custom_nstop`: Custom truncation point. Default is `0`, and the empirical formula 

```math
n_{\mathrm{stop}}=\max(x+4\sqrt[3]{x}+2, |m|x)
``` 

will be used.

Scattering information is outputed as a named tuple, including the following fields:

- `q_ext`: Extinction efficiency. Defined as ``C_{\mathrm{ext}}/\pi r^2``, where ``C_{\mathrm{ext}}`` is the extinction cross section.
- `q_sca`: Scattering efficiency.
- `q_abs`: Absorption efficiency.
- `q_back`: Backscattering efficiency. Defined as ``4(\mathrm{d}C_\mathrm{ext}/\mathrm{d}\Omega)/r^2``.
- `asymm`: Asymmetry factor ``\langle\cos(\theta)\rangle``.
- `S₁`, `S₂`: Amplitude scattering matrix components. See Eq. (3.12) in Bohren and Huffman (1983). Both S₁ and S₂ are vectors containing `Nθ` values.

References:

- Bohren, C.F., Huffman, D.R., 1983. Absorption and scattering of light by small particles. John Wiley & Sons.
"""
function bhmie(T, x::Real, m::Number, Nθ::Integer=181; nextra::Integer=15, custom_nstop::Integer=0)
    x = T(x)
    m = Complex{T}(m)
    y = m * x
    nstop = iszero(custom_nstop) ? Int(floor(Float64(max(x + 4 * ∛x + 2, x * abs(m))))) : custom_nstop
    nmax = nstop + nextra
    θ = collect(range(zero(T), T(π), Nθ))
    μ = cos.(θ)
    d = zeros(Complex{T}, nmax)
    for n in nmax-1:-1:1
        d[n] = (n + 1) / y - (1.0 / (d[n+1] + (n + 1) / y))
    end

    π_ = zeros(T, Nθ)
    τ = zeros(Nθ)
    π₀ = zeros(Nθ)
    π₁ = ones(Nθ)
    S₁ = zeros(Complex{T}, Nθ)
    S₂ = zeros(Complex{T}, Nθ)
    ψ₀ = cos(x)
    ψ₁ = sin(x)
    χ₀ = -sin(x)
    χ₁ = cos(x)
    ξ₁ = complex(ψ₁, -χ₁)
    a = zero(T)
    b = zero(T)
    q_sca = zero(T)
    asymm = zero(T)
    for n in 1:nstop
        fn = T((2n + 1) / (n * (n + 1)))
        ψ = (2n - 1) * ψ₁ / x - ψ₀
        χ = (2n - 1) * χ₁ / x - χ₀
        ξ = complex(ψ, -χ)
        aₙ = ((d[n] / m + n / x) * ψ - ψ₁)
        aₙ /= ((d[n] / m + n / x) * ξ - ξ₁)
        bₙ = ((d[n] * m + n / x) * ψ - ψ₁) / ((d[n] * m + n / x) * ξ - ξ₁)
        q_sca += (abs(aₙ)^2 + abs(bₙ)^2) * (2n + 1)
        asymm += (real(aₙ) * real(bₙ) + imag(aₙ) * imag(bₙ)) * fn
        if n > 1
            asymm += (real(a) * real(aₙ) + imag(a) * imag(aₙ) + real(b) * real(bₙ) + imag(b) * imag(bₙ)) * (n - 1) * (n + 1) / n
        end

        for i in 1:(Nθ+1)÷2
            ii = Nθ + 1 - i
            π_[i] = π₁[i]
            τ[i] = n * μ[i] * π_[i] - (n + 1) * π₀[i]
            S₁[i] += fn * (aₙ * π_[i] + bₙ * τ[i])
            S₂[i] += fn * (aₙ * τ[i] + bₙ * π_[i])
            if i < ii
                p = n & 1 == 0 ? 1 : -1
                S₁[ii] += fn * (aₙ * π_[i] - bₙ * τ[i]) * p
                S₂[ii] += fn * (-aₙ * τ[i] + bₙ * π_[i]) * p
            end
        end

        ψ₀, ψ₁ = ψ₁, ψ
        χ₀, χ₁ = χ₁, χ
        ξ₁ = complex(ψ₁, -χ₁)
        a, b = aₙ, bₙ

        @. π₁ = ((2n + 1) * μ * π_ - (n + 1) * π₀) / n
        π₀ .= π_
    end

    asymm *= 2 / q_sca
    q_sca *= 2 / x^2
    q_ext = 4 / x^2 * real(S₁[1])
    q_abs = q_ext - q_sca
    q_back = 4 / x^2 * abs(S₁[end])^2

    return (q_ext=q_ext, q_sca=q_sca, q_abs=q_abs, q_back=q_back, asymm=asymm, S₁=S₁, S₂=S₂)
end

bhmie(x::Real, m::Number, Nθ::Integer=181; kwargs...) = bhmie(Float64, x, m, Nθ; kwargs...)
