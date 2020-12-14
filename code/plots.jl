
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