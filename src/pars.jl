
using Parameters

@with_kw immutable Param @deftype Float64

    ρ = 0.05
    γ = 2.0
    invγ = 1./γ
    r = 0.03
    Δhjb = 1000.
end
