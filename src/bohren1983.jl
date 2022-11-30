export bhmie, bhcoat

@doc raw"""
`bhmie([T=Float64,], x, m; Nθ=181, nextra=15, custom_nstop=0)`

Inputs:

- `T`: Type used for calculation. All real numbers will be stored as `T`, while all complex numbers will be stored as `Complex{T}`.
- `x`: Size parameter of the sphere scatterer. Defined as ``\frac{2\pi r}{\lambda}``
- `m`: Relative refractive index of the scatterer.

Keyword arguments:

- `Nθ`: Number of scattering angles to calculate. Default is `181`.
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
- `S` = (`S₁`, `S₂`): Amplitude scattering matrix components. See Eq. (3.12) in Bohren and Huffman (1983). Both S₁ and S₂ are vectors containing `Nθ` values.
- `F` = (`F₁₁`, `F₁₂`, `F₃₃`, `F₃₄`): Mueller scattering matrix components. See Eq. (3.16) in Bohren and Huffman (1983). All Fᵢⱼ are vectors containing `Nθ` values.

References:

- Bohren, C.F., Huffman, D.R., 1983. Absorption and scattering of light by small particles. John Wiley & Sons.
"""
function bhmie(T, x, m; Nθ = 181, nextra = 15, custom_nstop = 0)
    x = T(x)
    m = Complex{T}(m)
    y = m * x
    nstop = iszero(custom_nstop) ? Int(floor(Float64(max(x + 4 * ∛x + 2, x * abs(m))))) :
            custom_nstop
    nmax = nstop + nextra
    θ = collect(range(zero(T), T(π), Nθ))
    μ = cos.(θ)
    d = zeros(Complex{T}, nmax)
    for n in (nmax - 1):-1:1
        d[n] = (n + 1) / y - (1.0 / (d[n + 1] + (n + 1) / y))
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
            asymm += (real(a) * real(aₙ) + imag(a) * imag(aₙ) + real(b) * real(bₙ) +
                      imag(b) * imag(bₙ)) * (n - 1) * (n + 1) / n
        end

        for i in 1:((Nθ + 1) ÷ 2)
            ii = Nθ + 1 - i
            π_[i] = π₁[i]
            τ[i] = n * μ[i] * π_[i] - (n + 1) * π₀[i]
            S₁[i] += fn * (aₙ * π_[i] + bₙ * τ[i])
            S₂[i] += fn * (aₙ * τ[i] + bₙ * π_[i])
            if i < ii
                p = n & 1 == 0 ? -1 : 1
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

    F₁₁ = zeros(T, Nθ)
    F₁₂ = zeros(T, Nθ)
    F₃₃ = zeros(T, Nθ)
    F₃₄ = zeros(T, Nθ)
    coeff = 4.0 / (q_sca * x^2)

    for i in 1:Nθ
        F₁₁[i] = (abs2(S₁[i]) + abs2(S₂[i])) * 0.5 * coeff
        F₁₂[i] = -(abs2(S₁[i]) - abs2(S₂[i])) * 0.5 * coeff
        F₃₃[i] = real(S₂[i] * conj(S₁[i])) * coeff
        F₃₄[i] = imag(S₂[i] * conj(S₁[i])) * coeff
    end

    return (q_ext = q_ext, q_sca = q_sca, q_abs = q_abs, q_back = q_back, asymm = asymm,
            S = (S₁, S₂), F = (F₁₁, F₁₂, F₃₃, F₃₄))
end

function bhmie(x, m; Nθ = 181, nextra = 15, custom_nstop = 0)
    bhmie(Float64, x, m; Nθ = Nθ, nextra = nextra, custom_nstop = custom_nstop)
end

@doc raw"""
`bhcoat([T=Float64,], x_core, x_mantle, m_core, m_mantle; Nθ = 181, tolerance = 1e-8)`

Inputs:

- `T`: Type used for calculation. All real numbers will be stored as `T`, while all complex numbers will be stored as `Complex{T}`.
- `x_core`: Size parameter of the core. Defined as ``\frac{2\pi r}{\lambda}``
- `x_mantle`: Size parameter of the coated sphere. `x_mantle >= x_core` should hold.
- `m_core`: Refractive index of the core, relative to the host medium.
- `m_mantle`: Refractive index of the mantle, relative to the host medium.

Keyword arguments:

- `Nθ`: Number of scattering angles to calculate. Default is `181`.
- `tolerance`: Error tolerance. Default is `1e-8`.

Scattering information is outputed as a named tuple, including the following fields:

- `q_ext`: Extinction efficiency. Defined as ``C_{\mathrm{ext}}/\pi r^2``, where ``C_{\mathrm{ext}}`` is the extinction cross section.
- `q_sca`: Scattering efficiency.
- `q_abs`: Absorption efficiency.
- `q_back`: Backscattering efficiency. Defined as ``4(\mathrm{d}C_\mathrm{ext}/\mathrm{d}\Omega)/r^2``.
- `asymm`: Asymmetry factor ``\langle\cos(\theta)\rangle``.
- `S` = (`S₁`, `S₂`): Amplitude scattering matrix components. See Eq. (3.12) in Bohren and Huffman (1983). Both S₁ and S₂ are vectors containing `Nθ` values.
- `F` = (`F₁₁`, `F₁₂`, `F₃₃`, `F₃₄`): Mueller scattering matrix components. See Eq. (3.16) in Bohren and Huffman (1983). All Fᵢⱼ are vectors containing `Nθ` values.

References:

- Bohren, C.F., Huffman, D.R., 1983. Absorption and scattering of light by small particles. John Wiley & Sons.
"""
function bhcoat(T, x_core, x_mantle, m_core, m_mantle; Nθ = 181, tolerance = 1e-8)
    @assert x_mantle>=x_core "x_mantle must be greater than or equal to x_core"

    θ = collect(range(zero(T), T(π), Nθ))
    μ = cos.(θ)
    π₀ = zeros(T, Nθ)
    π₁ = ones(T, Nθ)
    π_ = zeros(T, Nθ)
    τ = zeros(T, Nθ)
    S₁ = zeros(Complex{T}, Nθ)
    S₂ = zeros(Complex{T}, Nθ)
    x = T(x_core)
    y = T(x_mantle)
    m₁ = Complex{T}(m_core)
    m₂ = Complex{T}(m_mantle)
    x₁ = m₁ * x
    x₂ = m₂ * x
    y₂ = m₂ * y
    m = m₂ / m₁
    nstop = trunc(Int, y + 4.0 * ∛y + 2)
    a = zeros(Complex{T}, nstop)
    b = zeros(Complex{T}, nstop)
    an = zero(Complex{T})
    bn = zero(Complex{T})
    d0x1 = cot(x₁)
    d0x2 = cot(x₂)
    d0y2 = cot(y₂)
    d1x1 = zero(T)
    d1x2 = zero(T)
    brack = zero(T)
    crack = zero(T)
    psi0y = cos(y)
    psi1y = sin(y)
    chi0y = -sin(y)
    chi1y = cos(y)
    xi0y = psi0y - 1.0im * chi0y
    xi1y = psi1y - 1.0im * chi1y
    chi0y2 = -sin(y₂)
    chi1y2 = cos(y₂)
    chipy2 = zero(T)
    chiy2 = zero(T)
    chi0x2 = -sin(x₂)
    chi1x2 = cos(x₂)
    chipx2 = zero(T)
    chix2 = zero(T)
    qsca = zero(T)
    qext = zero(T)
    asymm = zero(T)
    xback = zero(Complex{T})
    amess1 = zero(T)
    amess2 = zero(T)
    amess3 = zero(T)
    amess4 = zero(T)
    iflag = false
    for n in 1:nstop
        rn = T(n)
        fn = (2rn + 1) / (rn * (rn + 1))
        en = rn
        psiy = (2rn - 1) * psi1y / y - psi0y
        chiy = (2rn - 1) * chi1y / y - chi0y
        xiy = psiy - 1.0im * chiy
        d1y2 = 1 / (rn / y₂ - d0y2) - rn / y₂
        if !iflag
            d1x1 = 1 / (rn / x₁ - d0x1) - rn / x₁
            d1x2 = 1 / (rn / x₂ - d0x2) - rn / x₂
            chix2 = (2rn - 1) * chi1x2 / x₂ - chi0x2
            chiy2 = (2rn - 1) * chi1y2 / y₂ - chi0y2
            chipx2 = chi1x2 - rn * chix2 / x₂
            chipy2 = chi1y2 - rn * chiy2 / y₂
            ancap = m * d1x1 - d1x2
            ancap = ancap / (m * d1x1 * chix2 - chipx2)
            ancap = ancap / (chix2 * d1x2 - chipx2)
            brack = ancap * (chiy2 * d1y2 - chipy2)
            bncap = m * d1x2 - d1x1
            bncap = bncap / (m * chipx2 - d1x1 * chix2)
            bncap = bncap / (chix2 * d1x2 - chipx2)
            crack = bncap * (chiy2 * d1y2 - chipy2)
            amess1 = brack * chipy2
            amess2 = brack * chiy2
            amess3 = crack * chipy2
            amess4 = crack * chiy2
        end
        if abs(amess1) < tolerance * abs(d1y2) && abs(amess2) < tolerance &&
           abs(amess3) < tolerance * abs(d1y2) && abs(amess4) < tolerance
            brack = zero(T)
            crack = zero(T)
            iflag = true
        else
            iflag = false
        end
        dnbar = d1y2 - brack * chipy2
        dnbar = dnbar / (1.0 - brack * chiy2)
        gnbar = d1y2 - crack * chipy2
        gnbar = gnbar / (1.0 - crack * chiy2)
        if n > 1
            an1 = an
            bn1 = bn
            a[n - 1] = -an
            b[n - 1] = -bn
        end
        an = ((dnbar / m₂ + rn / y) * psiy - psi1y) / ((dnbar / m₂ + rn / y) * xiy - xi1y)
        bn = ((m₂ * gnbar + rn / y) * psiy - psi1y) / ((m₂ * gnbar + rn / y) * xiy - xi1y)
        qsca += (2rn + 1) * (abs(an)^2 + abs(bn)^2)
        xback += (2rn + 1) * (-1)^n * (an - bn)
        qext += (2rn + 1) * (real(an) + real(bn))
        asymm += fn * (real(an) * real(bn) + imag(an) * imag(bn))
        if n > 1
            asymm += (en - 1) * (en + 1) / en * (real(an1) * real(an) +
                      imag(an1) * imag(an) +
                      real(bn1) * real(bn) + imag(bn1) * imag(bn))
        end
        for j1 in 1:((Nθ + 1) ÷ 2)
            jj = Nθ + 1 - j1
            π_[j1] = π₁[j1]
            τ[j1] = rn * μ[j1] * π_[j1] - (rn + 1) * π₀[j1]
            S₁[j1] += fn * (an * π_[j1] + bn * τ[j1])
            S₂[j1] += fn * (an * τ[j1] + bn * π_[j1])
            if j1 ≠ jj
                p = n & 1 == 0 ? -1 : 1
                S₁[jj] += fn * (an * π_[j1] - bn * τ[j1]) * p
                S₂[jj] += fn * (-an * τ[j1] + bn * π_[j1]) * p
            end
        end
        psi0y = psi1y
        psi1y = psiy
        chi0y = chi1y
        chi1y = chiy
        xi1y = psi1y - 1.0im * chi1y
        chi0x2 = chi1x2
        chi1x2 = chix2
        chi0y2 = chi1y2
        chi1y2 = chiy2
        d0x1 = d1x1
        d0x2 = d1x2
        d0y2 = d1y2
        for j1 in 1:Nθ
            π₁[j1] = (2rn + 1) / rn * μ[j1] * π_[j1] - (rn + 1) * π₀[j1] / rn
            π₀[j1] = π_[j1]
        end
    end

    asymm *= 2 / qsca
    q_sca = 2qsca / y^2
    q_ext = 2qext / y^2
    q_abs = q_ext - q_sca
    q_back = abs(xback)^2 / y^2

    F₁₁ = zeros(T, Nθ)
    F₁₂ = zeros(T, Nθ)
    F₃₃ = zeros(T, Nθ)
    F₃₄ = zeros(T, Nθ)
    coeff = 4.0 / (q_sca * x^2)

    for i in 1:Nθ
        F₁₁[i] = (abs2(S₁[i]) + abs2(S₂[i])) * 0.5 * coeff
        F₁₂[i] = -(abs2(S₁[i]) - abs2(S₂[i])) * 0.5 * coeff
        F₃₃[i] = real(S₂[i] * conj(S₁[i])) * coeff
        F₃₄[i] = imag(S₂[i] * conj(S₁[i])) * coeff
    end

    return (q_ext = q_ext, q_sca = q_sca, q_abs = q_abs, q_back = q_back, asymm = asymm,
            S = (S₁, S₂), F = (F₁₁, F₁₂, F₃₃, F₃₄))
end

function bhcoat(x_core, x_mantle, m_core, m_mantle; Nθ = 181, tolerance = 1e-8)
    bhcoat(Float64, x_core, x_mantle, m_core, m_mantle; Nθ = Nθ, tolerance = tolerance)
end
