using Distributed
using FileIO, JLD2, Statistics, StatsBase
@everywhere using SharedArrays, ShallowWaters
@everywhere set_zero_subnormals(true)

nxs = round.(Int,10.0 .^ (1.5:0.1:3.7))      # increase from nx=32 to almost 2512
ndays = 3e7 * (1 ./ nxs).^2.6                # number of days to integrate for
experiments = ["Float16",
                "Float16, uncompensated",
                "Float16/32",
                "Float32",
                "Float64"]

compensateds = [true,false,false,false,false]
Ts = [Float16,Float16,Float16,Float32,Float64]
Tprogs = [Float16,Float16,Float32,Float32,Float64]

nsamples = 5
stats = ["median","min","max"]
timings = fill(0.0,length(experiments),length(nxs),length(stats))

println("BENCHMARKING SHALLOWWATERS WITH $(nworkers()) PROCESSES.\n")

for (iT,(T,Tprog,compensated)) in enumerate(zip(Ts,Tprogs,compensateds))
    println("\n-- $T")
    println("Compile on every worker:")

    # trigger compilation
    @sync begin
        @async for p in nworkers()
            @spawnat p run_model(T,Ndays=1;Tprog,compensated)
        end
    end

    println("\nBenchmark now:")

    for (inx,(nx,Ndays)) in enumerate(zip(nxs,ndays))
        println("nx=$nx, Ndays=$Ndays")
        B = SharedVector{Float64}(nsamples)

        ngridpoints = 1/2*nx^2

        if ngridpoints*sizeof(T) > 3e7 # serial for very large grids
            @sync @spawnat 2 for isample in 1:nsamples
                GC.gc()
                B[isample] = run_model(T,Tprog=Tprog,nx=nx,Ndays=Ndays,compensated=compensated,return_time=true)
            end
        else    # parallel otherwise
            @sync @distributed for isample in 1:nsamples
                GC.gc()
                B[isample] = run_model(T,Tprog=Tprog,nx=nx,Ndays=Ndays,compensated=compensated,return_time=true)
            end
        end

        for (istat,stats) in enumerate([median,minimum,maximum])
            timings[iT,inx,istat] = stats(B)
        end
    end
end

@save "benchmarking.jld2" timings experiments nxs stats
println("Benchmarking stored.")