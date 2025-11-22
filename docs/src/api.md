# API Reference

```@autodocs
Modules = [SpheroidalWaveFunctions]
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
