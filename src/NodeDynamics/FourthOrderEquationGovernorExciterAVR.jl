# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# @doc doc"""
# ```Julia
# FourthEq(H, P, D, Ω, E_f, T_d_dash ,T_q_dash ,X_q_dash ,X_d_dash,X_d, X_q)
# ```
#
# A node type that applies the 4th-order synchronous machine model
# with frequency/angle and voltage dynamics.
#
# Additionally to ``u``, it has the internal dynamic variables
# * ``\omega`` representing the frequency of the rotator relative to the grid frequency ``\Omega``, i.e. the real frequency ``\omega_r`` of the rotator is given as ``\omega_r = \Omega + \omega`` and
# * ``\theta`` representing the relative angle of the rotor with respect to the voltage angle ``\phi``.
#
# # Keyword Arguments
# - `H`: inertia
# - `P`: active (real) power output
# - `D`: damping coefficient
# - `Ω`: rated frequency of the power grid, often 50Hz
# - `T_d_dash`: time constant of d-axis
# - `T_q_dash`: time constant of q-axis
# - `X_d_dash`: transient reactance of d-axis
# - `X_q_dash`: transient reactance of q-axis
# - `X_d`: reactance of d-axis
# - `X_d`: reactance of q-axis
#
# # Mathematical Representation
# Using `FourthEq` for node ``a`` applies the equations
# ```math
#     u = -je_c e^{j\theta} = -j(e_d + je_q)e^{j\theta}\\
#     e_c= e_d + je_q = jue^{-j\theta}\\
#     i  = -ji'e^{j\theta} = -j(i_d+ j i_q )e^{j\theta} = Y^L \cdot u \\
#     i_c= i_d + ji_q = jie^{-j\theta}\\
#     p = \Re (i^* u)
# ```
# The fourth-order equations read (according to Sauer, p. 140, eqs. (6110)-(6114)) and p. 35 eqs(3.90)-(3.91)
# ```math
#     \frac{d\theta}{dt} = \omega \\
#      \frac{d\omega}{dt} = P-D\omega - p -(x'_q-x'_d)i_d i_q\\
#     \frac{d e_q}{dt} = \frac{1}{T'_d} (- e_q - (x_d - x'_d) i_{d}+ e_f) \\
#     \frac{d e_d}{dt} = \frac{1}{T'_q} (- e_d + (x_q - x'_q) i_{q})  \\
# ```
# With the PowerDynamics.jl \time{naming conventions} of $i$ and $u$ they read as
# ```math
#    \dot u = \frac{d}{dt}(-j e_c e^{j\theta})=-j(\dot e_d + j\dot e_q)e^{j\theta} + uj\omega
# ```
# """

@DynamicNode FourthOrderEqGovernorExciterAVR(H, P, D, Ω, T_d_dash ,T_q_dash ,X_q_dash ,X_d_dash,X_d, X_q, T_e, T_a, T_f, K_e, K_a, K_f, V_ref, R_d, T_sv, T_ch) <: OrdinaryNodeDynamics() begin # S_e_fd,
    @assert H > 0 "inertia (H) should be >0"
    @assert D >= 0 "damping (D) should be >=0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert T_q_dash > 0 "time constant of q-axis (T_q_dash) should be >0"
    @assert X_d_dash >= 0 "transient reactance of d-axis (X_d_dash) should be >=0"
    @assert X_q_dash >= 0 "transient reactance of q-axis (X_q_dash) should be >=0"
    @assert X_d >= 0 "reactance of d-axis (X_d_dash) should be >=0"
    @assert X_q >= 0 "reactance of q-axis (X_q_dash) should be >=0"
    # @assert T_e > 0
    # @assert T_a > 0
    # @assert T_f > 0
    # @assert K_e > -10
    # @assert K_a > 0
    # @assert K_f > 0
    # @assert V_ref > 0

    Ω_H = (Ω * 2pi) / H #    norm = 2 * H / (2 * np.pi * 50)  # normalize the parameters as done for coupling_const, input_power, damping_const

end [[θ, dθ],[ω, dω],[e_f, de_f],[v_r, dv_r],[r_f,dr_f],[P_sv, dP_sv],[P_m, dP_m]] begin
    ##### Model R1 #####
    i_c = 1im*i*exp(-1im*θ)
    e_c = 1im*u*exp(-1im*θ)
    p = real(u * conj(i))
    e_d = real(e_c)
    e_q = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    #Exciter
    V_mes = e_c - 1im*X_d_dash*i_c
    dr_f = (1 / T_f) * (- r_f + ((K_f/T_f) * e_f))
    dv_r = (1 / T_a) * (- v_r + (K_a * r_f) - ((K_a * K_f)/T_f)*e_f + K_a*(V_ref - abs(V_mes)))
    de_f = (1 / T_e) * ((- (K_e + (0.098*exp(0.55*e_f))) * e_f) + v_r) #S_e_fd, Sauer p70

    #Governor

    dP_sv = (1/T_sv) * (-P_sv + P - (1/R_d)*(((ω+(Ω*2pi))/(Ω*2pi))-1))
    dP_m  = (1/T_ch) * (-P_m  + P_sv)

    #Synchronous machine
    dθ = ω
    de_d = (1 / T_q_dash) * (- e_d + (X_q - X_q_dash) * i_q)
    de_q = (1 / T_d_dash) * (- e_q - (X_d - X_d_dash) * i_d + e_f)
    de_c = de_d + 1im*de_q
    du = -1im*de_c*exp(1im*θ)+ u*1im*ω
    dω = (P_m - D*ω - p - (X_q_dash - X_d_dash) * i_d * i_q) * Ω_H
end

export FourthOrderEqGovernorExciterAVR
