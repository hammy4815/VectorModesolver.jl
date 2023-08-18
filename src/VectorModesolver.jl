module VectorModesolver

import Arpack: eigs 
import SparseArrays: spzeros
import Interpolations: interpolate, Linear, Gridded, extrapolate, Flat

include("Modesolver.jl")
include("Mode.jl")

export VectorialModesolver
export assemble
export solve

end # module VectorModesolver
