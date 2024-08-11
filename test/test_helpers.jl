using Test
using KiteModels

@testset "Testing KiteModels...." begin
    path=pwd()
    tmpdir=mktempdir()
    println("tmpdir: ", tmpdir)
    mkpath(tmpdir)
    cd(tmpdir)
    KiteModels.copy_examples()
    println("Examples: ", readdir(joinpath(tmpdir, "examples")))
    @test isfile(joinpath(tmpdir, "examples", "bench.jl"))
    @test isfile(joinpath(tmpdir, "examples", "compare_kps3_kps4.jl"))
    @test isfile(joinpath(tmpdir, "examples", "menu.jl"))
    @test isfile(joinpath(tmpdir, "examples", "reel_out_1p.jl"))
    @test isfile(joinpath(tmpdir, "examples", "reel_out_4p.jl"))
    @test isfile(joinpath(tmpdir, "examples", "reel_out_4p_torque_control.jl"))
    @test isfile(joinpath(tmpdir, "examples", "simulate_simple.jl"))
    @test isfile(joinpath(tmpdir, "examples", "simulate_steering.jl"))
    if ! Sys.iswindows()
        rm(tmpdir, recursive=true)
    end
    cd(path)
end
