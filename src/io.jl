# Use data file to extract mesh filename from attributes
function extract_mesh_filename(filename::String)
  # Open file for reading
  h5open(filename, "r") do file
    # Extract filename relative to data file
    mesh_file = read(attributes(file)["mesh_file"])

    return joinpath(dirname(filename), mesh_file)
  end
end


# Read in mesh file and return relevant data
function read_meshfile(filename::String)
  # Open file for reading
  h5open(filename, "r") do file
    # Extract basic information
    if exists(attributes(file), "ndims")
      ndims = read(attributes(file)["ndims"])
    else
      ndims = read(attributes(file)["ndim"]) # FIXME once Trixi's 3D branch is merged & released
    end
    n_children_per_cell = 2^ndims
    n_cells = read(attributes(file)["n_cells"])
    n_leaf_cells = read(attributes(file)["n_leaf_cells"])
    center_level_0 = read(attributes(file)["center_level_0"])
    length_level_0 = read(attributes(file)["length_level_0"])

    # Extract coordinates, levels, child cells
    coordinates = Array{Float64}(undef, ndims, n_cells)
    coordinates .= read(file["coordinates"])
    levels = Array{Int}(undef, n_cells)
    levels .= read(file["levels"])
    child_ids = Array{Int}(undef, n_children_per_cell, n_cells)
    child_ids .= read(file["child_ids"])

    # Extract leaf cells (= cells to be plotted) and contract all other arrays accordingly
    leaf_cells = similar(levels)
    n_cells = 0
    for cell_id in 1:length(levels)
      if sum(child_ids[:, cell_id]) > 0
        continue
      end

      n_cells += 1
      leaf_cells[n_cells] = cell_id
    end
    leaf_cells = leaf_cells[1:n_cells]

    coordinates = coordinates[:, leaf_cells]
    levels = levels[leaf_cells]

    return center_level_0, length_level_0, leaf_cells, coordinates, levels, ndims
  end
end


function read_datafile(filename::String)
  # Open file for reading
  h5open(filename, "r") do file
    ndims = read(attributes(file)["ndims"])
    # Extract basic information
    if exists(attributes(file), "polydeg")
      polydeg = read(attributes(file)["polydeg"])
    else
      polydeg = read(attributes(file)["N"])
    end
    n_elements = read(attributes(file)["n_elements"])
    n_variables = read(attributes(file)["n_vars"])
    time = read(attributes(file)["time"])

    # Extract labels for legend
    labels = Array{String}(undef, 1, n_variables)
    for v = 1:n_variables
      labels[1, v] = read(attributes(file["variables_$v"])["name"])
    end

    # Extract data arrays
    n_nodes = polydeg + 1

    if ndims == 3
      # Read 3d data
      data = Array{Float64}(undef, n_nodes, n_nodes, n_nodes, n_elements, n_variables)
    elseif ndims == 2
      # Read 2d data
      data = Array{Float64}(undef, n_nodes, n_nodes, n_elements, n_variables)
    else
      error("unsupported number of dimensions: $ndims")
    end
    for v = 1:n_variables
      vardata = read(file["variables_$v"])
      @views data[.., :, v][:] .= vardata
    end

    return labels, data, n_nodes, time
  end
end
