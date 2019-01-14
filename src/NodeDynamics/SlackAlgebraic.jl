# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
SlackAlgebraic(;U)
```

A node type that locally fixes the complex voltage (``U``) of the node.

As the complex voltage can be represented as ``U=Ve^{i\phi}``, this is equivlant
to fixing the voltage magnitude ``V`` and the angle ``\phi``.

# Keyword Arguments
- `U`: the complex voltage

# Mathematical Representation
Using `SlackAlgebraic` for node ``a`` applies the equation
```math
0 = U_a - u_a.
```
"""
@DynamicNode SlackAlgebraic(U) <: OrdinaryNodeDynamicsWithMass(m_u=false, m_int=no_internal_masses) begin
end [] begin
        du = u - U
end

getStaticParameters(n::SlackAlgebraic) = AbstractStaticSlack(V=abs(n.U),φ=angle(n.U))

export SlackAlgebraic
