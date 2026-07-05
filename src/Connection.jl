# src/Connection.jl

abstract type AbstractConnection end

struct LeviCivitaConnection <: AbstractConnection
end

christoffel(
    conn::LeviCivitaConnection,
    metric::FisherMetric,
    ρ,
    basis
)
