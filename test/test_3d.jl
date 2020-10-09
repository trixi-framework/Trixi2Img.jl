module Test3D

using Test
using Trixi2Img

include("test_trixi2img.jl")

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)


@testset "3D" begin
  run_trixi(joinpath("3d", "parameters.toml"), n_steps_max=1)

  @testset "uniform mesh as PNG" begin
    @testset "default slice" begin
      test_trixi2img("solution_000000.h5", outdir)
        # To be able to compare hashes, we would need to make sure that exactly the same versions
        # of all required libraries are installed everywhere, which does not seem to be a good option
        # right now.
    end

    @testset "x-axis slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:x, slice_axis_intersect=0.6)
    end

    @testset "y-axis slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:y, slice_axis_intersect=-0.35)
    end

    @testset "z-axis slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:z, slice_axis_intersect=-0.89)
    end

    @testset "negative border slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:y, slice_axis_intersect=-1)
    end

    @testset "positive border slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:z, slice_axis_intersect=1)
    end

    @testset "outside negative border slice" begin
      @test_throws ErrorException("slice_axis_intersect outside of domain") trixi2img(
          joinpath(outdir, "solution_000000.h5"); output_directory=outdir,
          slice_axis=:x, slice_axis_intersect=-1.01)
    end

    @testset "outside positive border slice" begin
      @test_throws ErrorException("slice_axis_intersect outside of domain") trixi2img(
          joinpath(outdir, "solution_000000.h5"); output_directory=outdir,
          slice_axis=:y, slice_axis_intersect=1.005)
    end
  end

  # PDF tests do not work on Windows
  if !Sys.iswindows()
    @testset "uniform mesh as PDF with grid lines" begin
      test_trixi2img("solution_000000.h5", outdir,
          format=:pdf, grid_lines=true)
          # To be able to compare hashes, we would need to make sure that exactly the same versions
          # of all required libraries are installed everywhere, which does not seem to be a good option
          # right now.
    end
  end
end

# Clean up afterwards: delete Trixi output directory
@test_nowarn rm(outdir, recursive=true)

end
