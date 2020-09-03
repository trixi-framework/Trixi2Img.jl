"""
    convert(filename::AbstractString...;
            format=:png, variables=[], verbose=false, grid_lines=false,
            output_directory=".", nvisnodes=nothing, max_supported_level=11)

Convert two-dimensional Trixi-generated output files to image files (PNG or PDF).

# Arguments
- `filename`: One or more 2D solution/restart files generated by Trixi to create an image from.
              Filenames support file globbing, e.g., "solution*" to match all files starting
              with `solution`.
- `format`: Output format for solution/restart files. Can be 'vtu' or 'vti'.
- `variables`: Names of the variables to create images for. If empty, each variable found in the
               input file(s) will be plotted.
- `verbose`: Set to `true` to enable verbose output.
- `grid_lines`: Plot outline of elements. (warning: this is an expensive operation!)
- `output_directory`: Output directory where generated files are stored.
- `nvisnodes`: Number of visualization nodes per element (default: twice the number of DG nodes).
               A value of `0` (zero) uses the number of nodes in the DG elements.
- `max_supported_level`: Maximum cell refinement level supported for plotting.

# Examples
```julia
julia> Trixi2Img.convert("out/solution_000*.h5")
[...]
```
"""
function Trixi2Img.convert(filename::AbstractString...;
                           format=:png, variables=[], verbose=false, grid_lines=false,
                           output_directory=".", nvisnodes=nothing, max_supported_level=11)
  # Note: We have to prefix `Trixi2Img.` to `convert` to avoid overwriting `Base.convert`
  # Reset timer
  reset_timer!()

  # Convert filenames to a single list of strings
  if isempty(filename)
    error("no input file was provided")
  end
  filenames = String[]
  for pattern in filename
    append!(filenames, glob(pattern))
  end

  # Ensure valid format
  if !(format in (:png, :pdf))
    error("unsupported output format '$format' (must be 'vtu' or 'vti')")
  end

  # Iterate over input files
  for filename in filenames
    verbose && println("Processing file $filename...")

    # Check if data file exists
    if !isfile(filename)
      error("data file '$filename' does not exist")
    end

    # Get mesh file name
    meshfile = extract_mesh_filename(filename)

    # Check if mesh file exists
    if !isfile(meshfile)
      error("mesh file '$meshfile' does not exist")
    end

    # Read mesh
    verbose && println("| Reading mesh file...")
    @timeit "read mesh" (center_level_0, length_level_0,
                         leaf_cells, coordinates, levels) = read_meshfile(meshfile)

    # Read data
    verbose && println("| Reading data file...")
    @timeit "read data" labels, unstructured_data, n_nodes, time = read_datafile(filename)

    # Determine resolution for data interpolation
    max_level = maximum(levels)
    if max_level > max_supported_level
      error("Maximum refinement level in data file $max_level is higher than " *
            "maximum supported level $max_supported_level")
    end
    max_available_nodes_per_finest_element = 2^(max_supported_level - max_level)
    if nvisnodes == nothing
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

    # Normalize element coordinates: move center to (0, 0) and domain size to [-1, 1]²
    n_elements = length(levels)
    normalized_coordinates = similar(coordinates)
    for element_id in 1:n_elements
      @views normalized_coordinates[:, element_id] .= (
           (coordinates[:, element_id] .- center_level_0) ./ (length_level_0 / 2 ))
    end

    # Interpolate unstructured DG data to structured data
    verbose && println("| Interpolating data...")
    @timeit "interpolate data" (structured_data =
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

    # Set up plotting
    gr()
    if format == :pdf
      GR.inline("pdf")
    elseif format == :png
      GR.inline("png")
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

    # Create output directory if it does not exist
    mkpath(output_directory)

    # Calculate x and y limits
    xlims = (-1, 1) .* (length_level_0/2 + center_level_0[1])
    ylims = (-1, 1) .* (length_level_0/2 + center_level_0[2])

    verbose && println("| Creating plots...")
    @timeit "plot generation" for (variable_id, label) in enumerate(variables)
      verbose && println("| | Variable $label...")

      # Create plot
      verbose && println("| | | Creating figure...")
      @timeit "create figure" plot(size=(2000,2000), thickness_scaling=1,
                                   aspectratio=:equal, legend=:none, title="$label (t = $time)",
                                   colorbar=true, xlims=xlims, ylims=ylims,
                                   tickfontsize=18, titlefontsize=28)

      # Plot contours
      verbose && println("| | | Plotting contours...")
      @timeit "plot contours" contourf!(xs, ys, node_centered_data[:, :, variable_id],
                                        c=:bluesreds, levels=128, linewidth=0)

      # Plot grid lines
      if grid_lines
        verbose && println("| | | Plotting grid lines...")
        @timeit "plot grid lines" plot!(vertices_x, vertices_y, linecolor=:black,
                                        linewidth=1, grid=false)
      end

      # Determine output file name
      base, _ = splitext(splitdir(filename)[2])
      output_filename = joinpath(output_directory, "$(base)_$(label)." * string(format))

      # Save file
      @timeit "save plot" savefig(output_filename)
    end
  end

  verbose && println("| done.\n")
  print_timer()
  println()
end

