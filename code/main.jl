import Pkg
Pkg.activate("code")

begin
    import CSV
    using DataFrames
    using EcologicalNetworks
    using LinearAlgebra, Statistics
    using Plots
    using ProgressMeter
    using EcologicalNetworksPlots
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
predictions = filter(x -> x.strength > R[x.from, x.to], predictions)
plot(
    [x.strength for x in predictions],
    c=:grey, lab="",
    fill = (:grey, 0, 0.2),
    frame = :origin,
    xlab = "Rank", ylab="Score"
    )

# Plot of additional predictions
I = initial(RandomInitialLayout, N)
@showprogress for step in 1:1000
  position!(ForceDirectedLayout(0.005), I, N)
end
plot(I, N, aspectratio=1)
scatter!(I, N, bipartite=true, nodesize=degree(N))

core3 = collect(keys(filter(p -> p.second >= 3, degree(N))))
Ncore = simplify(N[core3∩species(N; dims=1), core3∩species(N; dims=2)])

plot(I, Ncore)
scatter!(I, Ncore, nodesize=degree(N), bipartite=true)

K = BipartiteNetwork(zeros(Bool, size(N)), EcologicalNetworks.species_objects(N)...)
for pre in predictions[1:100]
    K[pre.from, pre.to] = true
end
simplify!(K)

plot!(I, K, lc=:red)
scatter!(I, N, nodesize=degree(N), bipartite=true, mf=:white)
