using Test: @test_nowarn, @test
using SHA
using Trixi
using Trixi2Img

# pathof(Trixi) returns /path/to/Trixi/src/Trixi.jl, dirname gives the parent directory
const EXAMPLES_DIR = joinpath(pathof(Trixi) |> dirname |> dirname, "examples")


function run_trixi(elixir; parameters...)
  @test_nowarn trixi_include(joinpath(EXAMPLES_DIR, elixir); parameters...)
end


function sha1file(filename)
  open(filename) do f
    bytes2hex(sha1(f))
  end
end


function set_pdf_creation_date(filename, date="20200101080000")
  # Implementation based on https://discourse.julialang.org/t/removing-the-first-line-in-a-text-file/35521/12
  write(filename,
        open(filename) do f
          content = read(f, String)
          replace(content, r"/CreationDate \(D:[0-9]*\)"=>"/CreationDate (D:$date)");
        end)
end


function test_trixi2img(filenames, outdir; hashes=nothing, kwargs...)
  @test_nowarn trixi2img(joinpath(outdir, filenames); output_directory=outdir, kwargs...)

  if !isnothing(hashes)
    for (filename, hash_expected) in hashes
      # If output format is PDF, we need to reset the CreationDate to be able to compare hashes
      if endswith(filename, ".pdf")
        set_pdf_creation_date(joinpath(outdir, filename))
      end

      hash_measured = sha1file(joinpath(outdir, filename))
      @test hash_expected == hash_measured
    end
  end
end
