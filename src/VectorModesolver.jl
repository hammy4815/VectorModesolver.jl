module VectorModesolver

import Arpack: eigs 
import SparseArrays: spzeros, SparseMatrixCSC
import Interpolations: interpolate, Linear, Gridded, extrapolate, Flat
using Parameters
using Trapz
using Revise

export VectorialModesolver, assemble, solve, εtype

include("Modesolver.jl")
include("Mode.jl")
include("Visualization.jl")

end # module VectorModesolver
