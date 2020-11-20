module Test3D

using Test
using Trixi2Img

include("test_trixi2img.jl")

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)


@testset "3D" begin
  run_trixi(joinpath("3d", "elixir_advection_extended.jl"), maxiters=1)

  @testset "uniform mesh as PNG" begin
    @testset "default slice" begin
      test_trixi2img("solution_000000.h5", outdir)
        # To be able to compare hashes, we would need to make sure that exactly the same versions
        # of all required libraries are installed everywhere, which does not seem to be a good option
        # right now.
    end

    @testset "x-axis slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:x, slice_axis_intercept=0.6)
    end

    @testset "y-axis slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:y, slice_axis_intercept=-0.35)
    end

    @testset "z-axis slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:z, slice_axis_intercept=-0.89)
    end

    @testset "negative border slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:y, slice_axis_intercept=-1)
    end

    @testset "positive border slice" begin
      test_trixi2img("solution_000000.h5", outdir, slice_axis=:z, slice_axis_intercept=1)
    end

    @testset "outside negative border slice" begin
      slice_axis_intercept=-1.01
      slice_axis=:x
      @test_throws ErrorException(string(
          "Slice plane $slice_axis = $slice_axis_intercept outside of domain. ",
          "$slice_axis must be between -1.0 and 1.0")) trixi2img(
          joinpath(outdir, "solution_000000.h5"); output_directory=outdir,
          slice_axis=slice_axis, slice_axis_intercept=slice_axis_intercept)
    end

    @testset "outside positive border slice" begin
      slice_axis_intercept=1.005
      slice_axis=:y
      @test_throws ErrorException(string(
          "Slice plane $slice_axis = $slice_axis_intercept outside of domain. ",
          "$slice_axis must be between -1.0 and 1.0")) trixi2img(
          joinpath(outdir, "solution_000000.h5"); output_directory=outdir,
          slice_axis=slice_axis, slice_axis_intercept=slice_axis_intercept)
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
