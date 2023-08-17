module VectorModesolver

import Arpack: eigs 
import SparseArrays: spzeros
import Interpolations: interpolate, Linear, Gridded

include("Modesolver.jl")
include("Mode.jl")

end # module VectorModesolver
