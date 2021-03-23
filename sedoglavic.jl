using Symbolics, ModelingToolkit, LinearAlgebra, AbstractAlgebra

function linear_ode_local_identifiability(
    equations::Array,
    output_functions::Array,
    x_variables::Array{Num,1}, 
    input_functions::Array,
    parameters::Array,
    power_series_precision::Int=1
)
    @parameters t
    poly_ring, t_ring = PolynomialRing(QQ, string(t))
    n = length(x_variables)
    ℓ = length(parameters)

    

    u = sum([rand((1:100)) * t^(i - 1) for i in 0:(n + ℓ + 1)]) # random power series modulo t^(n+l+1)

    odes = [eqn.rhs - eqn.lhs for eqn in equations]

    # calculate jacobian wrt dx/dt for each x
    Jac_odes_wrt_x_dot = Symbolics.jacobian(odes, D.(x_variables))

    # calculate jacobian wrt each x
    Jac_odes_wrt_x = Symbolics.jacobian(odes, x_variables) 

    current_precision = 1

    while current_precision < power_series_precision
        # specialize parameters to random values
        # this, ideally, should be done
        


        current_precision *= 2
    end

end