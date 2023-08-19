module VectorModesolver

import Arpack: eigs 
import SparseArrays: spzeros, SparseMatrixCSC
import Interpolations: interpolate, Linear, Gridded, extrapolate, Flat
using Parameters
using Revise

export VectorialModesolver, assemble, solve, ϵtype

include("Modesolver.jl")
include("Mode.jl")
include("Visualization.jl")

end # module VectorModesolver
