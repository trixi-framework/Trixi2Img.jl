module TestManual

using Test
using Documenter
using Trixi2Img

DocMeta.setdocmeta!(Trixi2Img, :DocTestSetup, :(using Trixi2Img); recursive=true)
doctest(Trixi2Img, manual=false)

end # module
