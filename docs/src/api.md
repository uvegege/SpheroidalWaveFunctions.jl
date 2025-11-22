# API Reference

## Angular Functions

### Prolate Angular Functions

```@docs
prolate_angular_leg
prolate_angular_ps
```

### Oblate Angular Functions

```@docs
oblate_angular_leg
oblate_angular_ps
```

## Radial Functions

### Prolate Radial Functions

```@docs
prolate_radial1
prolate_radial2
```

### Oblate Radial Functions

```@docs
oblate_radial1
oblate_radial2
```

## Characteristic Values

### Single Eigenvalue

```@docs
prolate_cv
oblate_cv
```

### Eigenvalue Sequences

```@docs
prolate_cv_seq
oblate_cv_seq
```

## Notes

All angular and radial functions return tuples `(value, derivative)`:
- Angular functions: `(S, dS/dx)`
- Radial functions: `(R, dR/dξ)`

### Parameters

- `m::Integer`: Azimuthal order (m ≥ 0)
- `n::Integer`: Mode number (n ≥ m)
- `c::Number`: Spheroidal parameter
- `λ::Number`: Characteristic value (optional, computed if not provided)
- `dr::Vector`: Expansion coefficients (optional, computed if not provided)
- `c2k::Vector`: Power series coefficients (optional, computed if not provided)
- `x::Number`: Angular argument in [-1, 1]
- `ξ::Number`: Radial argument (ξ ≥ 1)

### Coordinate Systems

**Prolate spheroidal coordinates**: Used for elongated geometries.

**Oblate spheroidal coordinates**: Used for flattened geometries.
