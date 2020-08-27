module Trixi2Img

# Include other packages
using Glob: glob
using HDF5: h5open, attrs, exists
using Plots: plot, plot!, gr, savefig, contourf!
using TimerOutputs
import GR

# Number of spatial dimensions
const ndim = 2

# Maximum level of cells supported for plotting
const max_supported_level = 11 # -> at most 2^11 = 2048 visualization nodes

# Include all source files
include("interpolation.jl")
include("interpolate.jl")
include("io.jl")

# Include top-level convert method
include("convert.jl")

end # module Trixi2Img

