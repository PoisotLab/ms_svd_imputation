import CSV
using DataFrames
using EcologicalNetworks
using LinearAlgebra, Statistics, SparseArrays
using Plots
using ProgressMeter
using EcologicalNetworksPlots

# Read the data
const virionette = DataFrame(CSV.File("./data/virionette.csv"));

# Make a network
hosts, viruses = sort(unique(virionette.host_species)), sort(unique(virionette.virus_genus))
i = indexin(virionette.virus_genus, viruses)
j = indexin(virionette.host_species, hosts)
V = fill(true, length(i))
N = BipartiteNetwork{Bool,String}(sparse(i,j,V), viruses, hosts)

include("./src/functions.jl")

T = BipartiteQuantitativeNetwork(float.(copy(N.edges)), species(N; dims=1), species(N; dims=2))
O = copy(T)
R = EcologicalNetworks.linearfilter(N; α=[0.0, 1.0, 1.0, 1.0])

A = Array(N.edges)
nj = sortperm(vec(sum(A; dims=1)))
ni = sortperm(vec(sum(A; dims=2)))

plot(
    heatmap(Array(lowrank(N; r=1).edges)[ni, nj], c=:cividis, frame=:none, cbar=false),
    heatmap(Array(lowrank(N; r=3).edges)[ni, nj], c=:cividis, frame=:none, cbar=false),
    heatmap(Array(lowrank(N; r=10).edges)[ni, nj], c=:cividis, frame=:none, cbar=false),
    heatmap(Array(lowrank(N; r=60).edges)[ni, nj], c=:cividis, frame=:none, cbar=false),
    size=(1000, 1000),
    dpi = 300
)
savefig("../figures/lowrank_illustration.png")

@showprogress for i in eachindex(N.edges)
    if !(N.edges[i])
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

core3 = collect(keys(filter(p -> p.second >= 5, degree(N))))
Ncore = simplify(N[core3∩species(N; dims=1), core3∩species(N; dims=2)])

plot(I, Ncore)
scatter!(I, Ncore, nodesize=degree(N), bipartite=true)

K = BipartiteNetwork(zeros(Bool, size(N)), EcologicalNetworks.species_objects(N)...)
for pre in predictions[1:10]
    K[pre.from, pre.to] = true
end
simplify!(K)

plot!(I, K, lc=:red)
scatter!(I, N, nodesize=degree(N), bipartite=true, mf=:white)
