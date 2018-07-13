function output_ini(u,v,η)
    if output == 1
        xudim = NcDim("x",nux,values=x_u)
        yudim = NcDim("y",nuy,values=y_u)
        xvdim = NcDim("x",nvx,values=x_v)
        yvdim = NcDim("y",nvy,values=y_v)
        xTdim = NcDim("x",nx,values=x_T)
        yTdim = NcDim("y",ny,values=y_T)
        tdim = NcDim("t",nout_total,unlimited=true)
        #println("NC Dimensions created.")

        uvar = NcVar("u",[xudim,yudim,tdim],t=Float32)
        vvar = NcVar("v",[xvdim,yvdim,tdim],t=Float32)
        ηvar = NcVar("eta",[xTdim,yTdim,tdim],t=Float32)
        #println("NC variables created.")

        ncu = NetCDF.create(runpath*"u.nc",uvar,mode=NC_NETCDF4)
        ncv = NetCDF.create(runpath*"v.nc",vvar,mode=NC_NETCDF4)
        ncη = NetCDF.create(runpath*"eta.nc",ηvar,mode=NC_NETCDF4)
        #println("NC ncfiles created.")

        NetCDF.putvar(ncu,"u",Float32.(u),start=[1,1,1],count=[-1,-1,1])
        NetCDF.putvar(ncv,"v",Float32.(v),start=[1,1,1],count=[-1,-1,1])
        NetCDF.putvar(ncη,"eta",Float32.(η),start=[1,1,1],count=[-1,-1,1])
        #println("Initial conditions written to file.")

        iout = 2    # counter for output time steps

        # also output scripts
        scripts_output()

        return (ncu,ncv,ncη),iout
    else
        return nothing, nothing
    end
end

function output_nc(ncs,u,v,η,i,iout)
    if i % nout == 0    # output only every nout time steps
        if output == 1
            NetCDF.putvar(ncs[1],"u",Float32.(u),start=[1,1,iout],count=[-1,-1,1])
            NetCDF.putvar(ncs[2],"v",Float32.(v),start=[1,1,iout],count=[-1,-1,1])
            NetCDF.putvar(ncs[3],"eta",Float32.(η),start=[1,1,iout],count=[-1,-1,1])
            #println("Time step $iout written to file.")
            iout += 1
        end
    end

    return ncs,iout
end

function output_close(ncs,progrtxt)
    if output == 1
        for nc in ncs
            NetCDF.close(nc)
        end
        println("All data stored.")
        write(progrtxt,"All data stored.")
        close(progrtxt)
    end
end

function get_run_id_path()
        if output == 1
                runlist = filter(x->startswith(x,"run"),readdir(outpath))
                existing_runs = [parse(Int,id[4:end]) for id in runlist]
                if length(existing_runs) == 0           # if no runfolder exists yet
                        runpath = outpath*"run0000/"
                        mkdir(runpath)
                        return 0,runpath
                else
                        run_id = maximum(existing_runs)+1
                        runpath = outpath*"run"*@sprintf("%04d",run_id)*"/"
                        mkdir(runpath)
                        return run_id,runpath
                end
        else
                return 0,"no runpath"
        end
end

const run_id,runpath = get_run_id_path()

function scripts_output()
        if output == 1
                # copy all files in juls main folder
                mkdir(runpath*"scripts")
                for juliafile in filter(x->endswith(x,".jl"),readdir())
                        cp(juliafile,runpath*"scripts/"*juliafile)
                end

                # and also in the src folder
                mkdir(runpath*"scripts/src")
                for juliafile in filter(x->endswith(x,".jl"),readdir("src"))
                        cp("src/"*juliafile,runpath*"scripts/src/"*juliafile)
                end
        end
end
