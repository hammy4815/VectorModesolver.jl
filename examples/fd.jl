using Revise
using VectorModesolver
using PyCall
using SparseArrays
EMpy = pyimport("EMpy")
Numpy = pyimport("numpy")

function convert_to_julia_sparse(python_sparse)
    # Get the COO format data
    rows = Int64.(python_sparse.row .+ 1)  # Convert to 1-based index for Julia
    cols = Int64.(python_sparse.col .+ 1)  # Convert to 1-based index for Julia
    data = python_sparse.data
    
    # Convert to Julia SparseMatrixCSC
    julia_sparse = sparse(rows, cols, data)
    return julia_sparse
end

function getApy(λ,x,y,epsfunc)
    x = Numpy.array(x)
    y = Numpy.array(y)
    solver = EMpy.modesolvers.FD.VFDModeSolver(λ, x, y, epsfunc, ("0","0","0","0"))
    mat = solver.build_matrix()
    return convert_to_julia_sparse(mat.tocoo())
end

# Define the domain
# Should return a tuple with (ϵxx, ϵxy, ϵyx, ϵyy, ϵzz)
function ϵ(x, y)
    return (1.0, 2.77, 2.77, 2.0, 3.0)
end

function epsfunc(x_, y_)
    xx, yy = Numpy.meshgrid(x_, y_)
    eps = Numpy.zeros((size(xx)..., 5))
    eps[1:end,1:end,1] .= 1.0
    eps[1:end,1:end,2] .= 2.77
    eps[1:end,1:end,3] .= 2.77
    eps[1:end,1:end,4] .= 2.0
    eps[1:end,1:end,5] .= 3.0
    return eps
end

# Parameters
λ = 1.55
x = 0:0.1:2.5
y = 0:0.1:2.5
neigs = 1
tol = 1e-8
boundary = (0,0,0,0)

# Define the modesolver
solver = VectorialModesolver(λ,x,y,boundary,ϵ)

# Solve for the modes
A = real.(assemble(solver))

Apy = getApy(λ,x,y,epsfunc)

sum(abs.(A-Apy))