include("Math.jl")

using Main.Math

function test()  
    println(ϕ(6.0) - ϕ(5.0, pi / 2))
    println(ϕ(6.0) + ϕ(5.0, pi / 2))
end

test()


