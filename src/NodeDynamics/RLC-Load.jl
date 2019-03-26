# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

@doc doc"""
```Julia
RLCLoad()
```
A node type that represents the RLC load model.

# Keyword Arguments
- `R`:
- `L`:
- `C`:


# Mathematical Representation
```math
	\dfrac{du_C}{dt} = \frac{1}{C}i_L(t)\\
    \dfrac{di_L}{dt} = -\frac{R}{L} i_L(t)+\frac{1}{L} u(t)
```

TODO: reference

"""
@DynamicNode RLCLoad(R,L,C) <: OrdinaryNodeDynamics()  begin
    @assert R > 0 "Resistance should be >0"
    @assert L > 0 "Inductance should be >0"
    @assert C > 0 "Capacitance should be >0"

end [[u_C, du_C],[i_L,di_L]] begin
    i_L = i #?????
    du_C = 1/C *i_L
    di_L = R/L*i_L+1/L*u
    du = du_C #????
end

export RLCLoad
