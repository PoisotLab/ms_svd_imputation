import Pkg
Pkg.activate("code")

begin
    import CSV
    using DataFrames
    using EcologicalNetworks
    using LinearAlgebra, Statistics
    using Plots
end