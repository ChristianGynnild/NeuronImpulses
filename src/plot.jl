using NeuronImpulses
#using NeuronImpulses:linear_cable_explicit_euler, linear_cable_implicit_euler, linear_cable_crank_nicolson
using Plots

x = NeuronImpulses.x
V0 = NeuronImpulses.V0


result = linear_cable_analytical()
Plots.heatmap(result)
Plots.savefig("figures/linear-cable-analytical")

V0 = result[:, 1]


mkpath("figures")
Plots.plot(x,V0)
Plots.savefig("figures/initial-condition")

result = linear_cable_explicit_euler(V0)
boundary = 3.924230158980217e300
replace!(result, Inf=>boundary, -Inf=>-boundary, Nothing=>boundary)
result[result.>= boundary] .= boundary
result[result.<= -boundary] .= -boundary

Plots.heatmap(result)
Plots.savefig("figures/linear-cable-solution-explicit-euler")

result = linear_cable_implicit_euler(V0)
Plots.heatmap(result)
Plots.savefig("figures/linear-cable-solution-implicit-euler")


result = linear_cable_crank_nicolson(V0)
Plots.heatmap(result)
Plots.savefig("figures/linear-cable-crank-nicolson")





#for i in 0:5
#    Plots.plot(x, result[:,round(Int, 1 + (N-1)/6*i)])
#    number = i+1
#    Plots.savefig("figures/linear-cable-equation-euler-explicit$number")
#end
