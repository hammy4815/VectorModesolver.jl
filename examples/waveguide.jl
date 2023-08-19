using VectorModesolver
using GeometryPrimitives

function main()
    
    # TODO define geometry with GeometryPrimitives

    λ = 1.55
    x = [i for i in 0:0.03:2.5]
    y = [i for i in 0:0.05:2.5]
    neigs = 1
    tol = 1e-8
    boundary = (0,0,0,0)
    solver = VectorialModesolver(λ,x,y,boundary,ϵ)
    modes = solve(solver, neigs, tol)

    plot_mode_fields(modes[1])
end

main()

