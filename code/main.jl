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