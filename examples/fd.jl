using Revise
using VectorModesolver
using PyCall
using SparseArrays
using PyPlot
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

    if (0.75 < x < 1.75) && (0.95 < y < 1.55)
        return (4.0, 0.0, 0.0, 4.0, 4.0)
    end

    return (1.0, 0.0, 0.0, 1.0, 1.0)
end

function epsfunc(x_, y_)
    xx, yy = Numpy.meshgrid(x_, y_)
    eps = Numpy.zeros((size(xx)..., 5))
    eps[1:end,1:end,1] .= 1.0
    eps[1:end,1:end,2] .= 0
    eps[1:end,1:end,3] .= 0
    eps[1:end,1:end,4] .= 1.0
    eps[1:end,1:end,5] .= 1.0
    return eps
end

# Parameters
λ = 1.55
x = [i for i in 0:0.05:2.5]
y = [i for i in 0:0.05:2.5]
neigs = 1
tol = 1e-8
boundary = (0,0,0,0)

# Define the modesolver
solver = VectorialModesolver(λ,x,y,boundary,ϵ)

# # Solve for the modes
# A = assemble(solver)

# Apy = getApy(λ,x,y,epsfunc)

# # sum(abs.(A-Apy))

modes = solve(solver, neigs, tol)

# Generating some dummy data for fields (replace with your real data)
Ex = real.(modes[1].Ex)
Ey = real.(modes[1].Ey)
Ez = imag.(modes[1].Ez)
Hx = real.(modes[1].Hx)
Hy = real.(modes[1].Hy)
Hz = imag.(modes[1].Hz)
xx = x * ones(length(y))'
yy = ones(length(x)) * y'
eps = ((x,y)->ϵ(x,y)[1]).(xx, yy)

PyPlot.figure(figsize=(10, 10)) # Create a 3x2 layout

# Create the heatmaps
PyPlot.subplot(3, 3, 1)
PyPlot.imshow(Ex, cmap="RdBu")
PyPlot.title("Ex")

PyPlot.subplot(3, 3, 2)
PyPlot.imshow(Ey, cmap="RdBu")
PyPlot.title("Ey")

PyPlot.subplot(3, 3, 3)
PyPlot.imshow(Ez, cmap="RdBu")
PyPlot.title("Ez")

PyPlot.subplot(3, 3, 4)
PyPlot.imshow(Hx, cmap="RdBu")
PyPlot.title("Hx")

PyPlot.subplot(3, 3, 5)
PyPlot.imshow(Hy, cmap="RdBu")
PyPlot.title("Hy")

PyPlot.subplot(3, 3, 6)
PyPlot.imshow(Hz, cmap="RdBu")
PyPlot.title("Hz")

PyPlot.subplot(3, 3, 7)
PyPlot.imshow(eps, cmap="Greys", label="eps")
PyPlot.title("eps")


PyPlot.savefig("/Users/ianhammond/GitHub/VectorModesolver/examples/ims.png")