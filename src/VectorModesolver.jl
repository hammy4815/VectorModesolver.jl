module VectorModesolver

import Arpack: eigs
import Interpolations: Flat, Gridded, Linear, extrapolate, interpolate
import SparseArrays: SparseMatrixCSC, spzeros

using GeometryPrimitves
using Parameters
using Revise

export VectorialModesolver, assemble, solve, Material, Geometry

include("Modesolver.jl")
include("Mode.jl")
include("Visualization.jl")

end # module VectorModesolver
