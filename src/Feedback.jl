@with_kw mutable struct Feedback
    t0::Float64=time()
    nans_detected::Book=false
    progress_txt::IOStream
    output::Bool
    i::Int=0
    nt::Int
end


"""Returns a human readable string representing seconds in terms of days, hours, minutes or seconds."""
function readable_secs(secs::Real)
    days = Int(floor(secs/3600/24))
    hours = Int(floor((secs/3600) % 24))
    minutes = Int(floor((secs/60) % 60))
    seconds = Int(floor(secs%3600%60))
    secs1f = @sprintf "%.1fs" secs%3600%60
    secs2f = @sprintf "%.2fs" secs%3600%60

    if days > 0
        return "$(days)d, $(hours)h"
    elseif hours > 0
        return "$(hours)h, $(minutes)min"
    elseif minutes > 0
        return "$(minutes)min, $(seconds)s"
    elseif seconds > 10
        return secs1f
    else
        return secs2f
    end
end

"""Estimates the total time the model integration will take."""
function duration_estimate(i,t,nt,progrtxt,P::Parameter)
    time_per_step = (time()-t) / (i-nadvstep)
    time_total = Int(round(time_per_step*nt))
    time_to_go = Int(round(time_per_step*(nt-i)))

    s1 = "Model integration will take approximately "*readable_secs(time_total)*","
    s2 = "and is hopefully done on "*Dates.format(now() + Dates.Second(time_to_go),Dates.RFC1123Format)

    println(s1)     # print inline
    println(s2)
    if output == 1  # print in txt
        write(progrtxt,"\n"*s1*"\n")
        write(progrtxt,s2*"\n")
        flush(progrtxt)
    end
end

"""Returns a boolean whether the prognostic variables contains a NaN."""
function nan_detection!(Prog::PrognosticVars,feedback::Feedback)

    #TODO include a check for Posits, are posits <: AbstractFloat?
    #TODO include check for tracer by other means than nan? (semi-Lagrange is unconditionally stable...)

    @unpack u,v,η,sst = Prog

    n_nan = sum(isnan.(u)) + sum(isnan.(v)) + sum(isnan.(η)) + sum(isnan.(sst))
    if n_nan > 0
        feeback.nans_detected = true
    end
end

"""Initialises the progress txt file."""
function feedback_init(S::ModelSetup)
    @unpack output = P
    @unpack nt = G

    if output
        txt = open(runpath*"progress.txt","w")
        s = "Starting Juls run $run_id on "*Dates.format(now(),Dates.RFC1123Format)
        println(s)
        write(txt,s*"\n")
        write(txt,"Juls will integrate $(Ndays)days at a resolution of $(nx)x$(ny) with Δ=$(Δ/1e3)km\n")
        write(txt,"Initial conditions are ")
        if initial_cond == "rest"
            write(txt,"rest.\n")
        else
            write(txt,"last time step of run $init_run_id.\n")
        end
        write(txt,"Boundary conditions are $bc_x with lbc=$lbc.\n")
        write(txt,"Numtype is "*string(Numtype)*".\n")
        write(txt,"Time steps are (Lin,Diff,Advcor,Lagr,Output)\n")
        write(txt,"$dtint, $(dtint*nstep_diff), $(dtint*nstep_advcor), $dtadvint, $(output_dt*3600)\n")
        write(txt,"\nAll data will be stored in $runpath\n")
    else
        println("Starting Juls on "*Dates.format(now(),Dates.RFC1123Format)*" without output.")
        txt = nothing
    end

    return Feedback(progress_txt=txt,output=output,nt=nt)
end

"""Feedback function that calls duration estimate, nan_detection and progress."""
#function feedback(u,v,η,sst,i,t,nt,nans_detected,progrtxt,P::Parameter)
function feedback(Prog::PrognosticVars,feedback::Feedback)
    # if i == nadvstep # measure time after tracer advection executed once
    #     t = time()
    # elseif i == min(2*nadvstep,nadvstep+50)
    #     # after the tracer advection executed twice or at least 50 steps
    #     duration_estimate(i,t,nt,progrtxt)
    # end
    #
    # if !nans_detected
    #     if i % nout == 0    # only check for nans when output is produced
    #         nans_detected = nan_detection(u,v,η,sst)
    #         if nans_detected
    #             println(" NaNs detected at time step $i")
    #             if output == 1
    #                 write(progrtxt," NaNs detected at time step $i")
    #                 flush(progrtxt)
    #             end
    #         end
    #     end
    # end

    if i > 100      # show percentage only after duration is estimated
        progress(i,nt,progrtxt,P)
    end

    return t,nans_detected
end

"""Finalises the progress txt file."""
function feedback_end(feedback::Feedback)
    @unpack output,t0,progress_txt = feedback

    s = " Integration done in "*readable_secs(time()-t0)*"."
    println(s)
    if output
        write(progrtxt,"\n"*s[2:end]*"\n")  # close txt file with last output
        flush(progrtxt)
    end
end

"""Converts time step into percent for feedback."""
function progress(feedback::Feedback)

    @unpack i,nt,progress_txt,output = feedback

    if ((i+1)/nt*100 % 1) < (i/nt*100 % 1)  # update every 1 percent steps.
        percent = Int(round((i+1)/nt*100))
        print("\r\u1b[K")
        print("$percent%")
        if output && (percent % 5 == 0) # write out only every 5 percent step.
            write(progrtxt,"\n$percent%")
            flush(progrtxt)
        end
    end
end
