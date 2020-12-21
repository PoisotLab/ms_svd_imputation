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

parameter_pos = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
parameters_val = parameters[parameter_pos,:]
a2, a3, a4 = parameters_val.a2, parameters_val.a3, parameters_val.a4

T = BipartiteQuantitativeNetwork(float.(copy(N.edges)), species(N; dims=1), species(N; dims=2))
O = copy(T)
R = EcologicalNetworks.linearfilter(N; α=[0.0, a2, a3, a4])

for pred_rank in 1:10
    try
        Threads.@threads for i in 1:size(N,1)
            for j in 1:size(N, 2)
                t = copy(T)
                t[i,j] = float(R[i,j])
                init = t[i,j]
                if init != 1.0
                    Δ = 1.0
                    iter = 1
                    while Δ > 1e-2
                        iter += 1
                        A = lowrank(t; r=pred_rank)
                        Δ = abs(A[i,j] - t[i,j])
                        t[i,j] = A[i,j]
                        iter ≥ 50 && break
                    end
                    # We only update if the value increased
                    if t[i,j] > init
                        O[i,j] = t[i,j]
                    end
                end
            end
        end

        predictions = filter(i -> T[i.from, i.to] != 1.0, interactions(O))
        sort!(predictions, by=(x) -> x.strength, rev=true)

        out_df = DataFrame(from = String[], to = String[], strength = Float64[], initial = Float64[])
        for i in predictions
            push!(out_df, (i.from, i.to, i.strength, R[i.from, i.to]))
        end
        out_df.score = out_df.strength ./ out_df.initial
        filter!(r -> r.strength != 1.0, out_df)

        ispath("results") || mkdir("results")
        CSV.write(joinpath("results", "case_$(parameter_pos)_rank_$(pred_rank).csv"), out_df)
    catch
        continue
    end
end