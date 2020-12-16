using DataFrames
import CSV
using ProgressMeter
using StatsPlots

hyperparameters = DataFrame(CSV.File(joinpath("data", "hyperparameters.csv")))
virionette = DataFrame(CSV.File(joinpath("data", "virionette.csv")))
fresnel = DataFrame(CSV.File(joinpath("data", "fresnel.csv")))

# List of new bat betacov-hosts
new_bats = filter(r -> !ismissing(r.Source), fresnel).Sp
old_betacov_hosts = filter(r -> ismissing(r.Source), fresnel).Sp

# Read an imputation file
function read_imputation(case::Int64)
    cstring = string(case)
    cfile = joinpath("results", "case_$(cstring).csv")

    isfile(cfile) || return nothing

    imputation_file = DataFrame(CSV.File(cfile))
    imputation_file.score = log.(imputation_file.strength./imputation_file.initial)
    imputation_bat = filter(r -> r.to in fresnel.Sp, imputation_file)
    imputation_bat_betacov = filter(r -> r.from == "Betacoronavirus", imputation_bat)
    
    sort!(imputation_bat_betacov, :score, rev=true)

    matched = indexin(new_bats, imputation_bat_betacov.to)
    filter!(!isnothing, matched)

    top_5 = count(matched .<= 5)
    top_10 = count(matched .<= 10)
    highest = minimum(matched)/length(imputation_bat_betacov.to)
    lowest = maximum(matched)/length(imputation_bat_betacov.to)
    proportion = length(matched)/length(new_bats)
    
    return (case, top_5, top_10, highest, lowest, proportion)
end

result_df = DataFrame(case = Int64[], top5 = Int64[], top10 = Int64[], best = Float64[], worse = Float64[], proportion = Float64[])

@showprogress for i in 1:size(hyperparameters, 1)
    out = read_imputation(i)
    if !isnothing(out)
        push!(result_df, out)
    end
end

hyperparameters.case = 1:size(hyperparameters, 1)

res = rightjoin(hyperparameters, result_df; on=:case)
CSV.write("results.csv", res)

# Ternary plot values
res.x = 0.5 .* (2 .* res.a3 .+ res.a4)
res.y = 0.5 .* sqrt(3.0) .* res.a4

@df res scatter(:x, :y, grid=false, aspectratio=1, frame=:none, mc=:top10, legend=false)
