# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
PVAlgebraic(;P,V)
```

A node type that locally fixes the active power (``P``) and the voltage magnitude (``V``) of the node.

# Keyword Arguments
- `P`: the active (real) power output
- `V`: voltage magnitude

# Mathematical Representation
Using `PVAlgebraic` for node ``a`` applies the equations
```math
0 = P_a - \Re\left(u_a \cdot i_a^*\right), \\
0 = V_a - \left\|u_a\right\|.
```
"""
@DynamicNode PVAlgebraic(P, V) <: OrdinaryNodeDynamicsWithMass(m_u=false, m_int=no_internal_masses) begin
end [] begin
    v = abs(u)
    p = real(u*conj(i))
    du = (v-V) + im*(p-P)
end

getStaticApproximation(n::PVAlgebraic) = StaticPV(P=n.P, V=n.V)

getInternalSteadyStateApproximation(::PVAlgebraic, u, i) = InternalSteadyState{PVAlgebraic}()

export PVAlgebraic
