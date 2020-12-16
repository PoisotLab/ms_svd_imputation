import CSV
using DataFrames
using EcologicalNetworks
using LinearAlgebra, Statistics, SparseArrays

# Read the data
const virionette = DataFrame(CSV.File("./data/virionette.csv"));
const parameters = DataFrame(CSV.File("./hyperparameters.csv"))

# Make a network
hosts, viruses = sort(unique(virionette.host_species)), sort(unique(virionette.virus_genus))
i = indexin(virionette.virus_genus, viruses)
j = indexin(virionette.host_species, hosts)
V = fill(true, length(i))
N = BipartiteNetwork{Bool,String}(sparse(i,j,V), viruses, hosts)

include("./src/functions.jl")

parameter_pos = get(ENV, "SLURM_ARRAY_TASK_ID", 1)
parameters_val = parameters[parameter_pos,:]
pred_rank = parameters_val.r
a2, a3, a4 = parameters_val.a2, parameters_val.a3, parameters_val.a4

T = BipartiteQuantitativeNetwork(float.(copy(N.edges)), species(N; dims=1), species(N; dims=2))
O = copy(T)
R = EcologicalNetworks.linearfilter(N; α=[0.0, a2, a3, a4])

Threads.@threads for i in eachindex(N)
    t = copy(T)
    t[i] = float(R[i])
    init = t[i]
    Δ = 1.0
    iter = 1
    while Δ > 1e-2
        iter += 1
        A = lowrank(t; r=pred_rank)
        Δ = abs(A[i] - t[i])
        t[i] = A[i]
        iter ≥ 50 && break
    end
    # We only update if the value increased
    if t[i] > init
        O[i] = t[i]
    end
end

predictions = filter(i -> !(N[i.from, i.to]), interactions(O))
sort!(predictions, by=(x) -> x.strength, rev=true)

out_df = DataFrame(from = String[], to = String[], strength = Float64[], initial = FLoat64[])
for i in predictions
    push!(out_df, (i.from, i.to, i.strength, R[i.from, i.to]))
end
ispath("results") || mkdir("results")
CSV.write(joinpath("results", "case_$(parameter_pos).csv"), out_df)