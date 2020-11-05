import Pkg
Pkg.activate("code")

begin
    import CSV
    using DataFrames
    using EcologicalNetworks
    using LinearAlgebra, Statistics
    using Plots
    using ProgressMeter
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

function impute!(output::T1, tmp::T1, template::T2, position; r=2, maxiter=50, tolerance=1e-2) where {T1 <: QuantitativeNetwork, T2 <: ProbabilisticNetwork}
    orig = tmp.A[position]
    tmp.A[position] = template.A[position]
    Δ = 1.0
    iter = 1
    while Δ > tolerance
        iter += 1
        A = lowrank(tmp; r=r)
        Δ = abs(A[position] - tmp[position])
        tmp.A[position] = A.A[position]
        iter ≥ maxiter && break
    end
    output.A[position] = tmp.A[position] # STORE!
    tmp.A[position] = orig # RESET!
end

T = BipartiteQuantitativeNetwork(float.(copy(N.A)), EcologicalNetworks.species_objects(N)...)
O = copy(T)
R = EcologicalNetworks.linearfilter(N; α=[0.0, 1.0, 1.0, 1.0])

@showprogress for i in eachindex(N.A)
    if !(N.A[i])
        impute!(O, T, R, i; r=2)
    end
end

predictions = filter(i -> !(N[i.from, i.to]), interactions(O))
sort!(predictions, by=(x) -> x.strength, rev=true)
plot(
    [x.strength for x in predictions],
    c=:grey, lab="",
    fill = (:grey, 0, 0.2),
    frame = :origin,
    xlab = "Rank", ylab="Score"
    )
