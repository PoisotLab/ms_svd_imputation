function LinearAlgebra.svd(N::T) where {T <: AbstractEcologicalNetwork}
    return svd(Array(N.edges))
end

function LinearAlgebra.rank(N::T) where {T <: AbstractEcologicalNetwork}
    return rank(Array(N.edges))
end

function Base.getindex(N::T, I::C) where {T <: AbstractEcologicalNetwork, C <: CartesianIndex}
    return N[I.I...]
end

function Base.setindex!(N::T, V::A, I::C) where {T <: AbstractEcologicalNetwork, A <: Any, C <: CartesianIndex}
    N[I.I...] = V
end

function Base.eachindex(N::T) where {T <: AbstractEcologicalNetwork}
    return eachindex(N.edges)
end

function lowrank(N::T; r::Integer=1) where {T <: AbstractEcologicalNetwork}
    r >= rank(N) && throw(ArgumentError("r ($(r)) is larger than the rank ($(rank(N))) of the network"))
    factor = svd(N)
    factor.S[(r+1):end] .= zero(eltype(factor.S))
    _rt = T <: AbstractUnipartiteNetwork ? UnipartiteQuantitativeNetwork : BipartiteQuantitativeNetwork
    return _rt(factor.U * Diagonal(factor.S) * factor.Vt, EcologicalNetworks._species_objects(N)...)
end

#=
function impute!(output::T1, tmp::T1, template::T2, position; r=2, maxiter=50, tolerance=1e-2) where {T1 <: QuantitativeNetwork, T2 <: ProbabilisticNetwork}
    orig = tmp[position]
    tmp[position] = float(template[position])
    Δ = 1.0
    iter = 1
    while Δ > tolerance
        iter += 1
        A = lowrank(tmp; r=r)
        Δ = abs(A[position] - tmp[position])
        tmp[position] = A[position]
        iter ≥ maxiter && break
    end
    output[position] = tmp[position] # STORE!
    tmp[position] = orig # RESET!
end

function impute(tmp::T1, template::T2, position; r=2, maxiter=50, tolerance=1e-2) where {T1 <: QuantitativeNetwork, T2 <: ProbabilisticNetwork}
    output = copy(tmp)
    impute!(output, tmp, template, position; r=r, maxiter=maxiter, tolerance=1e-2)
    return output
end
=#