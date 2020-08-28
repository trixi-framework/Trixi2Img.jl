using Documenter
using Trixi2Img

# Define module-wide setup such that the respective modules are available in doctests
DocMeta.setdocmeta!(Trixi2Img,
                    :DocTestSetup,
                    :(push!(LOAD_PATH, ".."); using Trixi2Img);
                    recursive=true)

# Make documentation
makedocs(
    # Specify modules for which docstrings should be shown
    modules = [Trixi2Img],
    # Set sitename
    sitename="Trixi2Img",
    # Provide additional formatting options
    format = Documenter.HTML(
        # Disable pretty URLs during manual testing
        prettyurls = get(ENV, "CI", nothing) == "true",
        # Explicitly add favicon as asset
        assets = ["assets/favicon.ico"],
        # Set canonical URL to GitHub pages URL
        canonical = "https://trixi-framework.github.io/Trixi2Img.jl/stable"
    ),
    # Explicitly specify documentation structure
    pages = [
        "Home" => "index.md",
        "Reference" => [
            "Trixi2Img" => "reference/trixi2img.md",
        ],
        "License" => "license.md"
    ]
)

deploydocs(
    repo = "github.com/trixi-framework/Trixi2Img.jl",
)
