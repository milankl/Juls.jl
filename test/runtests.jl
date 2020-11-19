using ShallowWaters
using Test

@testset "No Forcing" begin
    Prog = RunModel(Ndays=1,Fx0=0)
    @test all(Prog.u .== 0.0f0)
    @test all(Prog.v .== 0.0f0)
    @test all(Prog.η .== 0.0f0)
end

@testset "Boundary Conditions" begin
    Prog = RunModel(Ndays=1,bc="periodic")
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m

    Prog = RunModel(Ndays=1,bc="nonperiodic")
    @test all(abs.(Prog.η) .< 10)

    Prog = RunModel(Ndays=1,α=0)     # free-slip
    @test all(abs.(Prog.η) .< 10)
end

@testset "Continuity Forcing" begin
    Prog = RunModel(Ndays=1,surface_relax=true)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunModel(Ndays=1,surface_forcing=true)
    @test all(abs.(Prog.η) .< 10)
    @test all(Prog.u .!= 0.0f0)
end

@testset "Mixed Precision" begin
    Prog = RunModel(Float32,Tprog=Float64,Ndays=1)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunModel(Float16,Tprog=Float32,Ndays=1)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)
end

@testset "SSPRK2 tests" begin
    Prog = RunModel(time_scheme="SSPRK2",RKs=2,cfl=0.7,nstep_advcor=1)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunModel(time_scheme="SSPRK2",RKs=3,cfl=1.2,nstep_advcor=1)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunModel(time_scheme="SSPRK2",RKs=4,cfl=1.4,nstep_advcor=1)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunModel(time_scheme="SSPRK2",RKs=2,cfl=0.7,nstep_advcor=0)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunModel(time_scheme="SSPRK2",RKs=3,cfl=1.2,nstep_advcor=0)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunModel(time_scheme="SSPRK2",RKs=4,cfl=1.2,nstep_advcor=0)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)
end

@testset "SSPRK3 tests" begin
    Prog = RunModel(time_scheme="SSPRK3",cfl=1.2,nstep_advcor=1)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunModel(time_scheme="SSPRK3",cfl=1.2,nstep_advcor=0)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)
end