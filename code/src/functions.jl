function LinearAlgebra.svd(N::T) where {T <: AbstractEcologicalNetwork}
    return svd(N.A)
end

function LinearAlgebra.rank(N::T) where {T <: AbstractEcologicalNetwork}
    return rank(N.A)
end

function lowrank(N::T; r::Integer=1) where {T <: AbstractEcologicalNetwork}
    r >= rank(N) && throw(ArgumentError("r ($(r)) is larger than the rank ($(rank(N))) of the network"))
    factor = svd(N)
    factor.S[(r+1):end] .= zero(eltype(factor.S))
    _rt = T <: AbstractUnipartiteNetwork ? UnipartiteQuantitativeNetwork : BipartiteQuantitativeNetwork
    return _rt(factor.U * Diagonal(factor.S) * factor.Vt, EcologicalNetworks.species_objects(N)...)
end
