@doc doc"""
```Julia
CSIVoltagePT1(;τ_{DT1},τ_{P1},τ_{Q1},τ_{P2},τ_{Q2},,K_P,K_Q,V_r,P_r,Q_r)
```

A node type that applies the frequency and voltage droop control to control the frequency and
voltage dynamics.

Additionally to ``u``, it has the internal dynamic variable ``\omega`` representing
the frequency of the rotator relative to the grid frequency ``\Omega``, i.e. the
real frequency ``\omega_r`` of the rotator is given as ``\omega_r = \Omega + \omega``.

# Keyword Arguments
- `τ_{DT1}`: time constant DT1 element
- `τ_{P1}`: time constant 1st active power PT1 element
- `τ_{Q1}`: time constant 1st reactive PT1 element
- `τ_{P2}`: time constant 2nd active power PT1 element
- `τ_{Q2}`: time constant 2nd reactive power PT1 element
- `K_P`: droop constant frequency droop
- `K_Q`: droop constant voltage droop
- `V_r`: reference/ desired voltage
- `P_r`: active (real) power desired value
- `Q_r`: reactive (imag) power desired value


# Mathematical Representation
Using `CSI` for node ``a`` applies the equations
```math
\tau_{DT1}\dot{\omega}_{int} = -\omega_{int}+\dot{\phi}\\
P_{droop}= K_P (\omega_{int}-\omega_{ref})+P_{r}\\
Q_{droop}= K_Q (v- V_{r})+Q_{r}
\tau_{P1}\dot{P}_{PT1}=-P_{PT1}+P_{droop}\\
\tau_{P2}\dot{P}_{CSI}=-P_{CSI}+P_{PT1}\\
\tau_{Q1}\dot{Q}_{PT1}=-Q_{PT1}+Q_{droop}\\
\tau_{Q2}\dot{Q}_{CSI}=-Q_{CSI}+Q_{PT1}
```
"""
@DynamicNode CSIminimal(τ_DT1,τ_P1,τ_P2,τ_Q1,τ_Q2,K_P,K_Q,V_r,P_r,Q_r) <: OrdinaryNodeDynamics() begin
    @assert τ_DT1 > 0 "time constant  DT1 element >0"
    @assert τ_P1 > 0 "time constant 1st active power PT1 element should be >0"
    @assert τ_Q1 > 0 "time constant 1st reactive power PT1 element should be >0"
    @assert τ_P2 > 0 "time constant 2nd active power PT1 element should be >0"
    @assert τ_Q2 > 0 "time constant 2nd reactive power PT1 element should be >0"
    @assert K_Q > 0 "reactive power droop constant should be >0"
    @assert K_P > 0 "active power droop constant reactive power measurement should be >0"

end [[x, dx],[P_PT1,dP_PT1],[Q_PT1,dQ_PT1]] begin
    p = real(u*conj(i))
    q = imag(u*conj(i))
    v = abs(u)
    φ = angle(u)
    # https://de.wikipedia.org/wiki/Zustandsraumdarstellung#%C3%9Cbertragungsfunktion
    dx = 1/τ_DT1*x + φ
    ω_int = 1/(τ_DT1^2*2π)*x+1/(τ_DT1*2π)*φ
    #dω_int = 1/τ_DT1*(-ω_int + ω)
    dP_PT1 = 1/τ_P1*(-P_PT1+ K_P*ω_int+P_r)
    dp = 1/τ_P2*(-p + P_PT1)
    dQ_PT1 = 1/τ_Q1*(-Q_PT1+ K_Q*(v- V_r)+Q_r)
    dq = 1/τ_Q2*(-q+Q_PT1)
    ds = dp + 1im*dq
    du = ds/conj(i)
end

export CSIminimal
