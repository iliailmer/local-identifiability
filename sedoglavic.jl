using Symbolics, ModelingToolkit, LinearAlgebra, AbstractAlgebra

function linear_ode_local_identifiability(
    equations::Array,
    output_functions::Array,
    x_variables::Array{Num,1}, 
    input_functions::Array,
    param::Array
)
    @parameters t
    solution = Dict()
    n = length(x_variables)
    ℓ = length(param)
    ẋ = D.(x_variables)
    @variables ẋ_sub[1:length(x)](t) # placeholder for Differential(x[i])
    power_series_precision = n + ℓ + 1
    # random special value for x variables
    initial_conditions = Dict([x_variables[i] => rand((1:100)) for i in 1:length(x_variables)])
    for U in input_functions
        # random power series modulo t^(n+l+1)
        solution[U] = sum([rand((1:100)) * t^(i - 1) for i in 1:(n + ℓ + 1)])
    end
    odes = [eqn.rhs - eqn.lhs for eqn in equations]
    # calculate jacobian wrt dx/dt for each x
    Jac_odes_wrt_x_dot = Symbolics.jacobian(odes, ẋ)
    # calculate jacobian wrt each x
    Jac_odes_wrt_x = Symbolics.jacobian(odes, x_variables) 

    current_precision = 1
    power_series_ring, t_power_series = PowerSeriesRing(QQ, power_series_precision, string(t))

    polynomial_ring_generators = vcat(string.(x_variables), string.(param))
    polynomial_ring, ring_vars = PolynomialRing(QQ, polynomial_ring_generators)

    for x in x_variables
        solution[x] = power_series_ring(initial_conditions[x])
    end
    for x_dot in ẋ
        solution[x_dot] = power_series_ring(0)
        set_precision!(solution[x_dot], 1)
    end
    deriv2sub = Dict([(x_variables[i],  ẋ_sub[i]) for i in 1:length(x_variables)])
    sub2deriv = Dict([(ẋ_sub[i], x_variables[i]) for i in 1:length(x_variables)]) 
    
    # error occurs here
    ode_sub = [substitute(diff_eq, deriv2sub) for diff_eq in odes]
    
    while current_precision < power_series_precision
        # set precision of the current power series
        for i in 1:length(x_variables)
            set_prec!(solution[x_variables[i]], 2 * current_precision)
            set_prec!(solution[ẋ[i]], 2 * current_precision)
        end
        # specialize parameters to values
        eval_point = [set_precision!(solution[v], 2 * current_precision) for v in gens(polynomial_ring)]        
        current_precision *= 2
    end

end