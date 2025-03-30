using NeuronImpulses
#using NeuronImpulses:linear_cable_explicit_euler, linear_cable_implicit_euler, linear_cable_crank_nicolson
using Plots
using Statistics

x = NeuronImpulses.x
#V0 = NeuronImpulses.V0

println("Plotting analytical linear cable equation solution")
result_analytical = linear_cable_analytical()
Plots.heatmap(result_analytical)
Plots.savefig("figures/linear-cable-analytical")

V0 = result_analytical[:, 1]


mkpath("figures")
Plots.plot(x,V0)
Plots.savefig("figures/initial-condition")

println("Plotting explicit euler linear cable equation solution")
result_explicit_euler = linear_cable_explicit_euler(V0)
boundary = 3.924230158980217e300
replace!(result_explicit_euler, Inf=>boundary, -Inf=>-boundary, Nothing=>boundary)
result_explicit_euler[result_explicit_euler.>= boundary] .= boundary
result_explicit_euler[result_explicit_euler.<= -boundary] .= -boundary

Plots.heatmap(result_explicit_euler)
Plots.savefig("figures/linear-cable-solution-explicit-euler")

println("Plotting implicit euler linear cable equation solution")
result_implicit_euler = linear_cable_implicit_euler(V0)
Plots.heatmap(result_implicit_euler)
Plots.savefig("figures/linear-cable-solution-implicit-euler")

println("Plotting crank nicolson linear cable equation solution")
result_crank_nicolson = linear_cable_crank_nicolson(V0)
Plots.heatmap(result_crank_nicolson)
Plots.savefig("figures/linear-cable-crank-nicolson")


errors = hcat(
    abs.(mean(result_explicit_euler-result_analytical, dims=1)'),
    abs.(mean(result_implicit_euler-result_analytical, dims=1)'),
    abs.(mean(result_crank_nicolson-result_analytical, dims=1)')
)

println("Plotting error for implicit euler linear cable equation solution")
Plots.plot(LinRange(0, NeuronImpulses.time, NeuronImpulses.iterations), errors, label=["euler explicit" "euler implicit" "crank nicolson"])
Plots.ylims!(0,0.1)
xlabel!("Time")
ylabel!("Average error")
title!("Error of numerical methods for non linear cable equation")
Plots.savefig("figures/error")




#for i in 0:5
#    Plots.plot(x, result[:,round(Int, 1 + (N-1)/6*i)])
#    number = i+1
#    Plots.savefig("figures/linear-cable-equation-euler-explicit$number")
#end
