"""
    spheroidal_ang_1(m, n, c, x)

Smn(x, γ²) = ∑₀  d₂ᵣ(c) Psₘ₊₂ᵣᵐ(η)

Smn(x, γ²) = ∑₀  d₂ᵣ₊₁(c) Psₘ₊₂ᵣ₊₁ᵐ(η)

"""
function spheroidal_ang_1(m, n, dr, x)
    #λ =  find_eigenvalue(m, n, c)
    #dr = compute_dr2_mix(m, n, c, λ)
    v = iseven(n-m) ? 0 : 1
    S = 0.0
    dS = 0.0
    plms = assoc_legendre_Pm(m, m+2*(length(dr)-1)+v, x)

    for k in eachindex(dr)
        coef = dr[k]
        Plm  = plms[1][2*k-1+v]
        dPlm = plms[2][2*k-1+v]
        S += coef * Plm
        dS += coef * dPlm
    end
    s = iseven(m) ? 1 : -1
    return s*S, s*dS
end

"""
    spheroidal_ang_2(m, n, c, x)

Smn(x, γ²) = (-1)ᵐ (1 - x²)^(m/2) ∑ₖ₌₀ c₂ₖ (1 - x²)ᵏ

Smn(x, γ²) = (-1)ᵐ x(1 - x²)^(m/2) ∑ₖ₌₀ c₂ₖ (1 - x²)ᵏ
"""
function spheroidal_ang_2(m, n, c2k, x)
    #λ =  find_eigenvalue(m, n, c)
    #dr = compute_dr2_mix(m, n, c, λ)
    #c2k = compute_c2k(m, n, c, λ, dr)
    S = 0.0
    dS = 0.0
    for r in eachindex(c2k)
        k = r-1
        coef = c2k[r]
        S += coef* (1-x^2)^(k)
        dS += coef * (-2*k*x) * (1-x^2)^(k-1)
    end
    if iseven(n-m)
        dfactor = -m*x*(1-x^2)^(m/2-1)
        factor = (1-x^2)^(m/2)
        dS = (-1)^m * (factor * dS + dfactor * S)
        S = (-1)^m * factor * S
    else
        dfactor = (1-x^2)^(m/2 - 1) * (1 - (m+1) * x^2)
        factor = x*(1-x^2)^(m/2)
        dS = (-1)^m * (factor * dS + dfactor * S)
        S =  (-1)^m * factor * S
    end
    s = iseven(m) ? 1 : -1
    return s*S, s*dS
end
