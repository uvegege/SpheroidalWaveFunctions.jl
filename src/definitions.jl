"""
    prolate_angular_leg(m, n, c, x)
    prolate_angular_leg(m, n, c, λ, x)
    prolate_angular_leg(m, n, c, λ, dr, x)

Compute the prolate angular spheroidal wave function and its derivative using the Legendre expansion.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter
- `λ::Number`: characteristic value (optional, computed if not provided)
- `dr::Vector`: expansion coefficients (optional, computed if not provided)
- `x::Number`: argument in [-1, 1]

# Returns
- `S`: value of the angular function Sₘₙ(c, x)
- `dS`: derivative dSₘₙ/dx

# Examples
```julia
S, dS = prolate_angular_leg(1, 2, 0.5, 0.3)
```
"""
function prolate_angular_leg(m, n, c, x)
    λ =  find_eigenvalue(m, n, c)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_ang_1(m, n, dr, x)
end

function prolate_angular_leg(m, n, c, λ, x)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_ang_1(m, n, dr, x)
end

function prolate_angular_leg(m, n, c, λ, dr, x)
    return spheroidal_ang_1(m, n, dr, x)
end

"""
    oblate_angular_leg(m, n, c, x)
    oblate_angular_leg(m, n, c, λ, x)
    oblate_angular_leg(m, n, c, λ, c2k, x)

Compute the oblate angular spheroidal wave function and its derivative using the Legendre expansion.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter
- `λ::Number`: characteristic value (optional, computed if not provided)
- `c2k::Vector`: expansion coefficients (optional, computed if not provided)
- `x::Number`: argument in [-1, 1]

# Returns
- `S`: value of the angular function Sₘₙ(-ic, x)
- `dS`: derivative dSₘₙ/dx

# Examples
```julia
S, dS = oblate_angular_leg(1, 2, 0.5, 0.3)
```
"""
oblate_angular_leg(m, n, c, x) = prolate_angular_leg(m, n, im*c, x)
oblate_angular_leg(m, n, c, λ, x) = prolate_angular_leg(m, n, im*c, λ, x)
oblate_angular_leg(m, n, c, λ, c2k, x) = prolate_angular_leg(m, n, im*c, λ, c2k, x)


"""
    prolate_angular_ps(m, n, c, x)
    prolate_angular_ps(m, n, c, λ, x)
    prolate_angular_ps(m, n, c, λ, c2k, x)

Compute the prolate angular spheroidal wave function and its derivative using the power series expansion.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter
- `λ::Number`: characteristic value (optional, computed if not provided)
- `c2k::Vector`: expansion coefficients (optional, computed if not provided)
- `x::Number`: argument in [-1, 1]

# Returns
- `S`: value of the angular function Sₘₙ(c, x)
- `dS`: derivative dSₘₙ/dx

# Examples
```julia
S, dS = prolate_angular_ps(1, 2, 0.5, 0.3)
```
"""
function prolate_angular_ps(m, n, c, x)
    λ =  find_eigenvalue(m, n, c)
    dr = compute_dr2_mix(m, n, c, λ)
    c2k = compute_c2k(m, n, dr)
    return spheroidal_ang_2(m, n, c2k, x)
end

function prolate_angular_ps(m, n, c, λ, x)
    dr = compute_dr2_mix(m, n, c, λ)
    c2k = compute_c2k(m, n, dr)
    return spheroidal_ang_2(m, n, c2k, x)
end

function prolate_angular_ps(m, n, c, λ, c2k, x)
    return spheroidal_ang_2(m, n, c2k, x)
end

"""
    oblate_angular_ps(m, n, c, x)
    oblate_angular_ps(m, n, c, λ, x)
    oblate_angular_ps(m, n, c, λ, c2k, x)

Compute the oblate angular spheroidal wave function and its derivative using the power series expansion.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter
- `λ::Number`: characteristic value (optional, computed if not provided)
- `c2k::Vector`: expansion coefficients (optional, computed if not provided)
- `x::Number`: argument in [-1, 1]

# Returns
- `S`: value of the angular function Sₘₙ(-ic, x)
- `dS`: derivative dSₘₙ/dx

# Examples
```julia
S, dS = oblate_angular_ps(1, 2, 0.5, 0.3)
```
"""
oblate_angular_ps(m, n, c, x) = prolate_angular_ps(m, n, im*c, x)
oblate_angular_ps(m, n, c, λ, x) = prolate_angular_ps(m, n, im*c, λ, x)
oblate_angular_ps(m, n, c, λ, c2k, x) = prolate_angular_ps(m, n, im*c, λ, c2k, x)


"""
    prolate_radial1(m, n, c, ξ)
    prolate_radial1(m, n, c, λ, ξ)
    prolate_radial1(m, n, c, λ, dr, ξ)

Compute the prolate radial spheroidal wave function of the first kind and its derivative.

The radial function is regular at the origin and corresponds to the spherical Bessel function expansion.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter
- `λ::Number`: characteristic value (optional, computed if not provided)
- `dr::Vector`: expansion coefficients (optional, computed if not provided)
- `ξ::Number`: radial argument (ξ ≥ 1)

# Returns
- `R`: value of the radial function Rₘₙ⁽¹⁾(c, ξ)
- `dR`: derivative dRₘₙ⁽¹⁾/dξ

# Examples
```julia
R, dR = prolate_radial1(1, 2, 0.5, 1.5)
```
"""
function prolate_radial1(m, n, c, ξ)
    λ =  find_eigenvalue(m, n, c)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_1(m, n, c, dr, ξ)
end

function prolate_radial1(m, n, c, λ, ξ)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_1(m, n, c, dr, ξ)
end

function prolate_radial1(m, n, c, λ, dr, ξ)
    return spheroidal_rad_1(m, n, c, dr, ξ)
end

"""
    prolate_radial2(m, n, c, ξ)
    prolate_radial2(m, n, c, λ, ξ)
    prolate_radial2(m, n, c, λ, dr, ξ)

Compute the prolate radial spheroidal wave function of the second kind and its derivative.

The radial function is irregular at the origin and corresponds to the spherical Neumann function expansion.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter
- `λ::Number`: characteristic value (optional, computed if not provided)
- `dr::Vector`: expansion coefficients (optional, computed if not provided)
- `ξ::Number`: radial argument (ξ ≥ 1)

# Returns
- `R`: value of the radial function Rₘₙ⁽²⁾(c, ξ)
- `dR`: derivative dRₘₙ⁽²⁾/dξ

# Examples
```julia
R, dR = prolate_radial2(1, 2, 0.5, 1.5)
```
"""
function prolate_radial2(m, n, c, ξ)
    λ =  find_eigenvalue(m, n, c)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_2(m, n, c, dr, ξ)
end

function prolate_radial2(m, n, c, λ, ξ)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_2(m, n, c, dr, ξ)
end

function prolate_radial2(m, n, c, λ, dr, ξ)
    return spheroidal_rad_2(m, n, c, dr, ξ)
end



"""
    oblate_radial1(m, n, c, ξ)
    oblate_radial1(m, n, c, λ, ξ)
    oblate_radial1(m, n, c, λ, dr, ξ)

Compute the oblate radial spheroidal wave function of the first kind and its derivative.

The radial function is regular at the origin and corresponds to the modified spherical Bessel function expansion.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter
- `λ::Number`: characteristic value (optional, computed if not provided)
- `dr::Vector`: expansion coefficients (optional, computed if not provided)
- `ξ::Number`: radial argument (ξ ≥ 1)

# Returns
- `R`: value of the radial function Rₘₙ⁽¹⁾(-ic, iξ)
- `dR`: derivative dRₘₙ⁽¹⁾/dξ

# Examples
```julia
R, dR = oblate_radial1(1, 2, 0.5, 1.5)
```
"""
function oblate_radial1(m, n, c, ξ)
    λ =  find_eigenvalue(m, n, im*c)
    dr = compute_dr2_mix(m, n, im*c, λ)
    return spheroidal_rad_1(m, n, im*c, dr, im*ξ)
end

function oblate_radial1(m, n, c, λ, ξ)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_1(m, n, im*c, dr, im*ξ)
end

function oblate_radial1(m, n, c, λ, dr, ξ)
    return spheroidal_rad_1(m, n, im*c, dr, im*ξ)
end

"""
    oblate_radial2(m, n, c, ξ)
    oblate_radial2(m, n, c, λ, ξ)
    oblate_radial2(m, n, c, λ, dr, ξ)

Compute the oblate radial spheroidal wave function of the second kind and its derivative.

The radial function is irregular at the origin and corresponds to the modified spherical Neumann function expansion.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter
- `λ::Number`: characteristic value (optional, computed if not provided)
- `dr::Vector`: expansion coefficients (optional, computed if not provided)
- `ξ::Number`: radial argument (ξ ≥ 1)

# Returns
- `R`: value of the radial function Rₘₙ⁽²⁾(-ic, iξ)
- `dR`: derivative dRₘₙ⁽²⁾/dξ

# Examples
```julia
R, dR = oblate_radial2(1, 2, 0.5, 1.5)
```
"""
function oblate_radial2(m, n, c, ξ)
    λ =  find_eigenvalue(m, n, im*c)
    dr = compute_dr2_mix(m, n, im*c, λ)
    return spheroidal_rad_2(m, n, im*c, dr, im*ξ)
end

function oblate_radial2(m, n, c, λ, ξ)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_2(m, n, im*c, dr, im*ξ)
end

function oblate_radial2(m, n, c, λ, dr, ξ)
    return spheroidal_rad_2(m, n, im*c, dr, im*ξ)
end


"""
    prolate_cv(m, n, c)

Compute the characteristic value (eigenvalue) for prolate spheroidal wave functions.

The characteristic value λₘₙ(c) is the eigenvalue of the differential equation for the
angular spheroidal wave function.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter

# Returns
- `λ`: characteristic value λₘₙ(c)

# Examples
```julia
λ = prolate_cv(1, 2, 0.5)
```
"""
function prolate_cv(m, n, c)
    return find_eigenvalue(m, n, c)
end

"""
    oblate_cv(m, n, c)

Compute the characteristic value (eigenvalue) for oblate spheroidal wave functions.

The characteristic value λₘₙ(-ic) is the eigenvalue of the differential equation for the
angular spheroidal wave function in the oblate case.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter

# Returns
- `λ`: characteristic value λₘₙ(-ic)

# Examples
```julia
λ = oblate_cv(1, 2, 0.5)
```
"""
function oblate_cv(m, n, c)
    return find_eigenvalue(m, n, im*c)
end

"""
    prolate_cv_seq(m, n, c)

Compute a sequence of characteristic values for prolate spheroidal wave functions.

Returns the characteristic values λₘₘ(c), λₘ₊₁,ₘ(c), ..., λₙₘ(c) for all modes from m to n.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: maximum mode number (n ≥ m)
- `c::Number`: spheroidal parameter

# Returns
- Vector of characteristic values [λₘₘ(c), λₘ₊₁,ₘ(c), ..., λₙₘ(c)]

# Examples
```julia
λ_seq = prolate_cv_seq(1, 5, 0.5)  # Returns [λ₁₁, λ₂₁, λ₃₁, λ₄₁, λ₅₁]
```
"""
function prolate_cv_seq(m, n, c)
    return find_eigenvalue_seq(m, n, c)
end

"""
    oblate_cv_seq(m, n, c)

Compute a sequence of characteristic values for oblate spheroidal wave functions.

Returns the characteristic values λₘₘ(-ic), λₘ₊₁,ₘ(-ic), ..., λₙₘ(-ic) for all modes from m to n.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: maximum mode number (n ≥ m)
- `c::Number`: spheroidal parameter

# Returns
- Vector of characteristic values [λₘₘ(-ic), λₘ₊₁,ₘ(-ic), ..., λₙₘ(-ic)]

# Examples
```julia
λ_seq = oblate_cv_seq(1, 5, 0.5)  # Returns [λ₁₁, λ₂₁, λ₃₁, λ₄₁, λ₅₁]
```
"""
function oblate_cv_seq(m, n, c)
    return find_eigenvalue_seq(m, n, im*c)
end
