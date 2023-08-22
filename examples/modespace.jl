using VectorModesolver

function (ε::εtype)(x::Float64, y::Float64) 

    if (0.5 < x < 2.0) && (0.90 < y < 1.5)
        return (4.0, 0.0, 0.0, 4.0, 4.0)
    end

    return (1.0, 0.0, 0.0, 1.0, 1.0)
end

function main()
    # Solve Modes
    ε = εtype()
    λ = 1.55
    x = [i for i in 0:0.03:2.5]
    y = [i for i in 0:0.05:2.5]
    neigs = 3
    tol = 1e-8
    boundary = (0,0,0,0)
    solver = VectorialModesolver(λ,x,y,boundary,ε)
    modes = solve(solver, neigs, tol)

    # Ensure Modes are orthonormal
    for i in 1:neigs
        for j in 1:neigs
            println("Inner Product of Mode $i and Mode $j: ", VectorModesolver.innerproduct(modes[i], modes[j]))
        end
    end

    # Take arbitrary mode combination
    true_coeffs = rand(neigs)
    combo = sum([true_coeffs[i] * modes[i] for i in 1:neigs])

    # Recover the coefficients
    recovered_coeffs = zeros(neigs)
    for i in 1:neigs
        coeffs = VectorModesolver.modecoeffs((combo.Ex, combo.Ey, combo.Ez), (combo.Hx, combo.Hy, combo.Hz), modes[i])
        recovered_coeffs[i] = coeffs[1]
    end

    # Print the results
    println("True Coefficients: ", true_coeffs)
    println("Recovered Coefficients: ", recovered_coeffs)
end