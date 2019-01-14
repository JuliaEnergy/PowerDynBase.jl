# (C) 2018 authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

abstract type AbstractStaticParameters end
abstract type AbstractStaticPQ <: AbstractStaticParameters end
abstract type AbstractStaticPV <: AbstractStaticParameters end
abstract type AbstractStaticSlack <: AbstractStaticParameters end
struct StaticPQ <: AbstractStaticPQ
    P
    Q
    StaticPQ(;P, Q) = new(P, Q)
end
struct StaticPQ <: AbstractStaticPV
    P
    V
    StaticPQ(;P, V) = new(P, V)
end
struct StaticSlack <:AbstractStaticSlack
    V
    φ
    StaticPQ(;V, φ) = new(V, φ)
end

"""add doc string"""
function getStaticParameters end

getStaticParameters(n::AbstractNodeDynamics) =  n |> parametersof |> getStaticParameters
