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

R = EcologicalNetworks.linearfilter(N; α=[0.0, 1.0, 1.0, 1.0])

svd(R.A)