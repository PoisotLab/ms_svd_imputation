import Pkg
Pkg.activate("code")

begin
    import CSV
    using DataFrames
    using EcologicalNetworks
    using LinearAlgebra, Statistics
    using Plots
end

# Read the data
const virionette = DataFrame(CSV.File("./code/data/virionette.csv"));

# Make a network
hosts, viruses = unique(virionette.host_species), unique(virionette.virus_genus)
A = zeros(Bool, (length(viruses), length(hosts)))
N = BipartiteNetwork(A, viruses, hosts)

for record in eachrow(virionette)
    N[record.virus_genus, record.host_species] = true
end

include("./code/src/functions.jl")

function impute!(output::T1, tmp::T1, data::T2, template::T3, position; r=2, maxiter=50, tolerance=1e-2) where {T1 <: QuantitativeNetwork, T2 <: AbstractEcologicalNetwork, T3 <: ProbabilisticNetwork}
    tmp.A[position] = template.A[position]
    Δ = 1.0
    iter = 1
    while Δ > tolerance
        iter += 1
        A = lowrank(tmp; r=r)
        Δ = abs(A[position] - N[position])
        tmp.A[position] = A.A[position]
        iter ≥ maxiter && break
    end
    output.A[position] = tmp.A[position]
    tmp.A[position] = data.A[position] # RESET!
end

T = BipartiteQuantitativeNetwork(float.(copy(N.A)), EcologicalNetworks.species_objects(N)...)
O = copy(T)
R = EcologicalNetworks.linearfilter(N; α=[0.0, 1.0, 1.0, 1.0])

@time for i in eachindex(R.A)
    if !(N.A[i])
        @info i
        @time impute!(O, T, N, R, i)
    end
end