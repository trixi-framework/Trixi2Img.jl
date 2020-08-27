module Trixi2Img

# Include other packages
using Glob: glob
using HDF5: h5open, attrs, exists
using Plots: plot, plot!, gr, savefig, contourf!
using TimerOutputs
import GR

# Number of spatial dimensions
"""
    ndims

Number of spatial dimensions (= 2).
"""
const ndim = 2

# Include all source files
include("interpolation.jl")
include("interpolate.jl")
include("io.jl")

# Include top-level convert method
include("convert.jl")

end # module Trixi2Img

