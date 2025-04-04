import Plots
using SparseArrays
using LinearAlgebra
using ProgressBars

use_progressbar = true

if (!use_progressbar)
    ProgressBar = (x -> x)
end


a, b = -2, 4
N::Int = 101*5
Δx = (b-a)/N
iterations = 20000*5
time = 0.1*3
Δt = time/iterations
e = MathConstants.e
μ = (a+b)/2
standard_deviation = 0.01

α = Δt/Δx^2
λ = 1
τ = 1

analytical_start_time = 0.001

V_star = -40*10^-3 # Volt

gamma = 5

g(V) = (100/(1 + e^gamma*(V_star - V))) + 1/5



gaussian(x) = e^-((x-μ)^2/(2*standard_deviation^2))


x = LinRange(a,b,N)
V0 = gaussian.(x)


export linear_cable_explicit_euler, linear_cable_implicit_euler, linear_cable_crank_nicolson
export linear_cable_analytical

function linear_cable_explicit_euler(V0::AbstractVector{T}) where {T}
    k1 = λ^2*α/τ
    k2 = ((-2*λ^2*α - Δt)/τ + 1)

    A = spdiagm(
        -1=> ones(N-1)*k1, 
         0=> ones(N)*k2,
         1=> ones(N-1)*k1, 
                 )

    A[1, 2] = 2*k1
    A[N, N-1] = 2*k1

    result = Array{T}(undef, N, iterations)
    value = V0

    for i in ProgressBar(1:iterations)
        result[:,i] = value
        value = A*value
    end
    return result
end

function linear_cable_implicit_euler(V0::AbstractVector{T}) where {T}
    k3 = -λ^2*α/τ
    k4 = ((2*λ^2*α + Δt)/τ + 1)

    B = spdiagm(
        -1=> ones(N-1)*k3, 
         0=> ones(N)*k4,
         1=> ones(N-1)*k3, 
                 )

    #print(B)

    B[1, 2] = 2*k3
    B[N, N-1] = 2*k3

    B_inv = inv(Matrix(B))

    result = Array{T}(undef, N, iterations)
    value = V0

    println("Dims:$(size(value'))")
    println("Dims:$(size(B))")

    

    for i in ProgressBar(1:iterations)
        result[:,i] = value
        #value =  B \ value
        value =  B_inv * value
    end
    return result
end

function linear_cable_crank_nicolson(V0::AbstractVector{T}) where {T}
    k1 = λ^2*α/τ
    k2 = ((-2*λ^2*α - Δt)/τ + 1)

    A = spdiagm(
        -1=> ones(N-1)*k1, 
         0=> ones(N)*k2,
         1=> ones(N-1)*k1, 
    )

    A[1, 2] = 2*k1
    A[N, N-1] = 2*k1


    k3 = -λ^2*α/τ
    k4 = ((2*λ^2*α + Δt)/τ + 1)

    B = spdiagm(
        -1=> ones(N-1)*k3, 
         0=> ones(N)*k4,
         1=> ones(N-1)*k3, 
    )

    B[1, 2] = 2*k3
    B[N, N-1] = 2*k3

    C = inv(Matrix(B+I))*(A+I)

    result = Array{T}(undef, N, iterations)
    value = V0

    for i in ProgressBar(1:iterations)
        result[:,i] = value
        value = C*value
    end
    return result
end


function linear_cable_analytical()
    V(x,t) = normalisation_term/sqrt(4*pi*(λ^2/τ)*t)*exp(-(x-x0)^2/(4*λ^2/τ*t) - t/τ)

    normalisation_term = 1
    x0 = (a+b)/2


    t = analytical_start_time
    result = Array{Float32}(undef, N, iterations)

    for i in 1:iterations
        result[:,i] = V.(x,t)
        t += Δt
    end
    return result
end

