using RecipesBase
using DiffEqBase
import .Trixi

struct TrixiPlot{RealT, SolutionT, SemidiscretizationT<:Trixi.AbstractSemidiscretization}
  u::SolutionT
  semi::SemidiscretizationT
  time::RealT
end

TrixiPlot(u, prob::ODEProblem, time=NaN) = TrixiPlot(Trixi.wrap_array(u, prob.p), prob.p, time)
TrixiPlot(u, sol::ODESolution, time=NaN) = TrixiPlot(Trixi.wrap_array(u, sol.prob.p), sol.prob.p, time)
TrixiPlot(sol::ODESolution) = TrixiPlot(Trixi.wrap_array(sol.u[end], sol.prob.p), sol.prob.p, sol.t[end])
TrixiPlot(u, ::Trixi.SemidiscretizationEulerGravity) = error("not implemented")

function Base.show(io::IO, tp::TrixiPlot)
  print(io, "TrixiPlot(...)")
end

@recipe function f(tp::TrixiPlot; grid_lines=false, max_supported_level=11, nvisnodes=nothing,
                                  slice_axis=:z, slice_axis_intercept=0)
  # Extract basic information from solution
  mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(tp.semi)

  center_level_0 = mesh.tree.center_level_0
  length_level_0 = mesh.tree.length_level_0
  leaf_cells = Trixi.leaf_cells(mesh.tree)
  coordinates = mesh.tree.coordinates[:, leaf_cells]
  levels = mesh.tree.levels[leaf_cells]
  ndims = Trixi.ndims(mesh)

  labels = Array{String}(undef, 1, Trixi.nvariables(equations))
  labels[1, :] .= Trixi.varnames(Trixi.cons2cons, equations)
  if ndims == 3
    # Read 3d data
    data = Array{Float64}(undef, Trixi.nnodes(solver), Trixi.nnodes(solver), Trixi.nnodes(solver),
                          Trixi.nelements(solver, cache), Trixi.nvariables(equations))
  elseif ndims == 2
    # Read 2d data
    data = Array{Float64}(undef, Trixi.nnodes(solver), Trixi.nnodes(solver),
                          Trixi.nelements(solver, cache), Trixi.nvariables(equations))
  else
    error("unsupported number of dimensions: $ndims")
  end
  for variable in Trixi.eachvariable(equations)
    @views data[.., :, variable] .= tp.u[variable, .., :]
  end
  n_nodes = Trixi.nnodes(solver)
  time = 0.0

  # Calculate x and y limits
  xlims = (-1, 1) .* (length_level_0/2 + center_level_0[1])
  ylims = (-1, 1) .* (length_level_0/2 + center_level_0[2])

  # Determine axis labels
  if ndims == 3
    # Extract plot labels
    if slice_axis === :x
      xlabel = "y"
      ylabel = "z"
    elseif slice_axis === :y
      xlabel = "x"
      ylabel = "z"
    elseif slice_axis === :z
      xlabel = "x"
      ylabel = "y"
    else
      error("illegal dimension '$slice_axis', supported dimensions are :x, :y, and :z")
    end
  elseif ndims == 2
    xlabel = "x"
    ylabel = "y"
  end

  # Create plot
  # plot(size=(2000,2000), thickness_scaling=1,
  #      aspectratio=:equal, legend=:none, title="$label (t = $time)",
  #      colorbar=true, xlims=xlims, ylims=ylims,
  #      xlabel=xlabel, ylabel=ylabel,
  #      labelfontsize=18, tickfontsize=18,
  #      titlefontsize=28)

  # # Plot contours
  # contourf!(xs, ys, node_centered_data[:, :, variable_id], c=:bluesreds, levels=128, linewidth=0,
  #           match_dimensions=true)

  # # Plot grid lines
  # if grid_lines
  #   plot!(vertices_x, vertices_y, linecolor=:black, linewidth=1, grid=false)
  # end

  # Create actual plot series and set attributes
  # size --> ((2000, 2000))
  # thickness_scaling --> 1
  aspect_ratio --> :equal
  # match_dimensions --> true
  # legend --> :none
  # colorbar --> true
  title --> labels[1,1]
  xlims --> xlims
  ylims --> ylims
  xguide --> xlabel
  yguide --> ylabel

  seriestype --> :contour
  fill --> true
  # seriescolor --> :bluesreds
  # linewidth --> 0

  get_contour_data(center_level_0, length_level_0, leaf_cells, coordinates, levels, ndims, labels,
                   data, n_nodes, time)
end


function get_contour_data(center_level_0, length_level_0, leaf_cells, coordinates, levels, ndims,
                          labels, unstructured_data, n_nodes, time)
  variables = []
  grid_lines = false
  nvisnodes = nothing
  max_supported_level = 11
  slice_axis = :z
  slice_axis_intercept = 0

  # Determine resolution for data interpolation
  max_level = maximum(levels)
  if max_level > max_supported_level
    error("Maximum refinement level in data file $max_level is higher than " *
          "maximum supported level $max_supported_level")
  end
  max_available_nodes_per_finest_element = 2^(max_supported_level - max_level)
  if nvisnodes === nothing
    max_nvisnodes = 2 * n_nodes
  elseif nvisnodes == 0
    max_nvisnodes = n_nodes
  else
    max_nvisnodes = nvisnodes
  end
  nvisnodes_at_max_level = min(max_available_nodes_per_finest_element, max_nvisnodes)
  resolution = nvisnodes_at_max_level * 2^max_level
  nvisnodes_per_level = [2^(max_level - level)*nvisnodes_at_max_level for level in 0:max_level]
  # nvisnodes_per_level is an array (accessed by "level + 1" to accommodate
  # level-0-cell) that contains the number of visualization nodes for any
  # refinement level to visualize on an equidistant grid

  if ndims == 3
    (unstructured_data, coordinates, levels,
        center_level_0) = unstructured_2d_to_3d(unstructured_data,
        coordinates, levels, length_level_0, center_level_0, slice_axis,
        slice_axis_intercept)
  end

  # Normalize element coordinates: move center to (0, 0) and domain size to [-1, 1]²
  n_elements = length(levels)
  normalized_coordinates = similar(coordinates)
  for element_id in 1:n_elements
    @views normalized_coordinates[:, element_id] .= (
          (coordinates[:, element_id] .- center_level_0) ./ (length_level_0 / 2 ))
  end

  # Interpolate unstructured DG data to structured data
  (structured_data =
      unstructured2structured(unstructured_data, normalized_coordinates,
                              levels, resolution, nvisnodes_per_level))

  # Interpolate cell-centered values to node-centered values
  node_centered_data = cell2node(structured_data)

  # Determine axis coordinates for contour plot
  xs = collect(range(-1, 1, length=resolution+1)) .* length_level_0/2 .+ center_level_0[1]
  ys = collect(range(-1, 1, length=resolution+1)) .* length_level_0/2 .+ center_level_0[2]

  # Determine element vertices to plot grid lines
  if grid_lines
    vertices_x, vertices_y = calc_vertices(coordinates, levels, length_level_0)
  end

  # Check that all variables exist in data file
  if isempty(variables)
    append!(variables, labels)
  else
    for var in variables
      if !(var in labels)
        error("variable '$var' does not exist in the data file $filename")
      end
    end
  end

  variable_id = 1
  label = variables[variable_id]
  

  # Create plot
  # plot(size=(2000,2000), thickness_scaling=1,
  #      aspectratio=:equal, legend=:none, title="$label (t = $time)",
  #      colorbar=true, xlims=xlims, ylims=ylims,
  #      xlabel=xlabel, ylabel=ylabel,
  #      labelfontsize=18, tickfontsize=18,
  #      titlefontsize=28)

  # # Plot contours
  # contourf!(xs, ys, node_centered_data[:, :, variable_id], c=:bluesreds, levels=128, linewidth=0,
  #           match_dimensions=true)

  # # Plot grid lines
  # if grid_lines
  #   plot!(vertices_x, vertices_y, linecolor=:black, linewidth=1, grid=false)
  # end
  return xs, ys, transpose(node_centered_data[:, :, variable_id])
end
