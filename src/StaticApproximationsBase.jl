# (C) 2018 authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

abstract type AbstractStaticApproximation end

abstract type AbstractStaticPQ <: AbstractStaticApproximation end
abstract type AbstractStaticPV <: AbstractStaticApproximation end
abstract type AbstractStaticSlack <: AbstractStaticApproximation end

struct StaticPQ <: AbstractStaticPQ
    P
    Q
    StaticPQ(;P, Q) = new(P, Q)
end
struct StaticPV <: AbstractStaticPV
    P
    V
    StaticPV(;P, V) = new(P, V)
end
struct StaticSlack <:AbstractStaticSlack
    V
    φ
    StaticSlack(;V, φ) = new(V, φ)
end

# function Base.:(==)(s1::AbstractInternalState{N}, s2::AbstractInternalState{N}) where N
#     s1.vals == s2.vals
# end
function Base.:≈(s1::T, s2::T) where {T <: AbstractStaticApproximation}
    for f in fieldnames(T)
        if !(getfield(s1, f) ≈ getfield(s2, f))
            return false
        end
    end
    true
end

"""add doc string"""
function getStaticApproximation end

getStaticApproximation(n::AbstractNodeDynamics) =  n |> parametersof |> getStaticParameters

"""add doc string"""
function getInternalSteadyStateApproximation end

getInternalSteadyStateApproximation(n::AbstractNodeDynamics, u, i) =
    getInternalSteadyStateApproximation(parametersof(n), u, i)
