using DataFrames
import CSV

rnks = collect(1:8)
a2 = Float64[]
a3 = Float64[]
a4 = Float64[]

n = 15

for i in LinRange(0.0, 1.0, n)
    for j in LinRange(0.0, 1.0, n)
        for k in LinRange(0.0, 1.0, n)
            if i+j+k == 1.0
                push!(a2, i)
                push!(a3, j)
                push!(a4, k)
            end
        end
    end
end

df = DataFrame(r=Int64[], a2=Float64[], a3=Float64[], a4=Float64[])
for r in rnks
    for i in eachindex(a2)
        push!(df, (r, a2[i], a3[i], a4[i]))
    end
end

CSV.write("code/hyperparameters.csv", df)