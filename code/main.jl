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
