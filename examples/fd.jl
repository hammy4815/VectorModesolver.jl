using Revise
using VectorModesolver
using PyCall
using SparseArrays
using PyPlot
using LinearAlgebra
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
# Should return a tuple with (εxx, εxy, εyx, εyy, εzz)
function (ε::εtype)(x::Float64, y::Float64) 

    if (0.75 < x < 1.75) && (0.95 < y < 1.55)
        return (4.0, 0.0, 0.0, 4.0, 4.0)
    end

    return (1.0, 0.0, 0.0, 1.0, 1.0)
end
ε = εtype()

function epsfunc(x_, y_)
    eps = Numpy.zeros((length(x_), length(y_), 5))
    for (i, x) in enumerate(x_)
        for (j, y) in enumerate(y_)
            eps[i,j,1] = ε(x,y)[1]
            eps[i,j,2] = ε(x,y)[2]
            eps[i,j,3] = ε(x,y)[3]
            eps[i,j,4] = ε(x,y)[4]
            eps[i,j,5] = ε(x,y)[5]
        end
    end
    return eps
end

function getempymodes(λ,x,y,epsfunc,neigs, tol)
    x = Numpy.array(x)
    y = Numpy.array(y)
    return EMpy.modesolvers.FD.VFDModeSolver(λ, x, y, epsfunc, ("0","0","0","0")).solve(
        neigs, tol
    ).modes
end

# 1 < - - - - yj- - - - - >
# ^ [1. 1. 1. 1. 1. 1. 1.]
# | [1. 1. 1. 1. 1. 1. 1.]
# | [1. 1. 4. 4. 1. 1. 1.]
# | [1. 1. 4. 4. 1. 1. 1.]
# | [1. 1. 4. 4. 1. 1. 1.]
# xi[1. 1. 4. 4. 1. 1. 1.]
# | [1. 1. 1. 1. 1. 1. 1.]
# | [1. 1. 1. 1. 1. 1. 1.]
# | [1. 1. 1. 1. 1. 1. 1.]
# v [1. 1. 1. 1. 1. 1. 1.]
# 
#    W
# S  +  N
#    E

# Parameters
λ = 1.55
x = [i for i in 0:0.03:2.5]
y = [i for i in 0:0.05:2.5]
neigs = 1
tol = 1e-8
boundary = (0,0,0,0)

# Define the modesolver
solver = VectorialModesolver(λ,x,y,boundary,ε)

# # Solve for the modes
# A = assemble(solver)
# Apy = getApy(λ,x,y,epsfunc)
# w = (A - Apy) .> 1e-5
# @show findnz(w)[1]
# @show findnz(w)[2]
# sum(abs.(A-Apy))
# sum(abs.(A-Apy))

# open("warntype_output.txt", "w") do f
#     redirect_stdout(f) do
#         @code_warntype optimize=true assemble(solver)
#     end
# end
# a = assemble(solver)
# @code_warntype optimize=true ε(1.0,1.0)

modes = solve(solver, neigs, tol)
# emodes = getempymodes(λ,x,y,epsfunc,neigs, tol)
# # modes = solve(Apy, solver, neigs, tol)

# Generating some dummy data for fields (replace with your real data)
Ex = real.(modes[1].Ex)
Ey = real.(modes[1].Ey)
Ez = imag.(modes[1].Ez)
Hx = real.(modes[1].Hx)
Hy = real.(modes[1].Hy)
Hz = imag.(modes[1].Hz)
xx = x * ones(length(y))'
yy = ones(length(x)) * y'
eps = ((x,y)->ε(x,y)[1]).(xx, yy)

PyPlot.figure(figsize=(16, 6)) # Create a 3x2 layout

# Create the heatmaps
PyPlot.subplot(2, 4, 1)
PyPlot.imshow(Ex, cmap="RdBu")
PyPlot.title("Ex")
PyPlot.colorbar()

PyPlot.subplot(2, 4, 2)
PyPlot.imshow(Ey, cmap="RdBu")
PyPlot.title("Ey")
PyPlot.colorbar()

PyPlot.subplot(2, 4, 3)
PyPlot.imshow(Ez, cmap="RdBu")
PyPlot.title("Ez")
PyPlot.colorbar()

PyPlot.subplot(2, 4, 5)
PyPlot.imshow(Hx, cmap="RdBu")
PyPlot.title("Hx")
PyPlot.colorbar()

PyPlot.subplot(2, 4, 6)
PyPlot.imshow(Hy, cmap="RdBu")
PyPlot.title("Hy")
PyPlot.colorbar()

PyPlot.subplot(2, 4, 7)
PyPlot.imshow(Hz, cmap="RdBu")
PyPlot.title("Hz")
PyPlot.colorbar()

PyPlot.subplot(2, 4, 4)
PyPlot.imshow(eps', cmap="Greys", label="eps")
PyPlot.title("eps")


PyPlot.savefig("/Users/ianhammond/GitHub/VectorModesolver/examples/ims.png")