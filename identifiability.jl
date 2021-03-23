using ModelingToolkit, Symbolics, LinearAlgebra

succeed = true
@parameters t, θ[1:4]
@variables x[1:2](t), y, U(t)
D = Differential(t)

n = length(x)
ℓ = length(θ)
ν = 1
equations = [D(x[1]) ~ x[1]^2 * θ[1] + θ[2] * x[1] * x[2] + U, D(x[2]) ~ θ[3] * x[1]^2 + θ[4] * x[1] * x[2]]
outputs = [y ~ x[1]]

# constructing P 
@variables ẋ[1:length(x)](t)

subs = [substitute(expr.lhs, Dict([D(x[i]) => ẋ[i] for i in 1:length(x)])) ~ expr.rhs for expr in ode]

F = [x.rhs for x in ode]
G = [x.rhs for x in outputs]

P = [x.lhs - x.rhs for x in subs]

# symbolic part: computing up to ν derivatives of G wrt t
# then compute jacobian wrt x and θ

function LieDerivative(f)
    
end

dict_deriv2rhs = Dict([D(x[i]) => F[i] for i in 1:length(x)])

# constructing ∇P
dPdẋ = Symbolics.jacobian(P, ẋ)
dPdx = Symbolics.jacobian(P, x)
dPdθ = Symbolics.jacobian(P, θ)

Γinit = Matrix{Int64}(I, n, n)
Λinit = Matrix{Int64}(0 * I, n, ℓ)
Φinit = [rand((1:100)) for i in 1:6] # random integers
Θinit = [rand((1:100)) for i in 1:6] # random integers
u = sum([rand((1:100)) * t^i for i in 0:6]) # random power series

# How to generate power series solution?

