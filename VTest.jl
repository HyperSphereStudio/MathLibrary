include("Math.jl")

using Main.Math
using Main.Math.VPlot

pts = []

for i in 1.0:500.0
    push!(pts, ϕ(i) * ϕ(i, pi))
end

Main.Math.VPlot.plotv(pts)
