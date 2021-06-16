"""Coefficients for strong stability-preserving Runge-Kutta 3rd order.
From: KETCHESON, LOĆZI, AND PARSANI, 2014. INTERNAL ERROR PROPAGATION IN EXPLICIT RUNGE–KUTTA METHODS, 
SIAM J NUMER ANAL 52/5. DOI:10.1137/130936245"""
struct SSPRK3coeff{T<:AbstractFloat}
    n::Int
    s::Int
    kn::Int
    mn::Int
    Δt_Δn::T
    kna::T
    knb::T
    Δt_Δnc::T
end

"""Generator function for a SSPRK3coeff struct."""
function SSPRK3coeff{T}(P::Parameter,Δt_Δ::T) where T
    n = P.RKn
    s = n^2
    kn = n*(n+1) ÷ 2 + 1
    mn = (n-1)*(n-2) ÷ 2 + 1
    Δt_Δn = convert(T,Δt_Δ/(n^2-n))
    kna = convert(T,(n-1)/(2n-1))
    knb = convert(T,n/(2n-1))
    Δt_Δnc = convert(T,Δt_Δ/(n*(2n-1)))

    return SSPRK3coeff{T}(n,s,kn,mn,Δt_Δn,kna,knb,Δt_Δnc)
end

struct Constants{T<:AbstractFloat,Tprog<:AbstractFloat}

    # RUNGE-KUTTA COEFFICIENTS 2nd/3rd/4th order including timestep Δt
    RKaΔt::Array{Tprog,1}
    RKbΔt::Array{Tprog,1}
    Δt_Δs::Tprog            # Δt/(s-1) wher s the number of stages
    Δt_Δ::Tprog             # Δt/Δ - timestep divided by grid spacing
    Δt_Δ_half::Tprog        # 1/2 * Δt/Δ
    SSPRK3c::SSPRK3coeff    # struct containing all coefficients for SSPRK3

    # BOUNDARY CONDITIONS
    one_minus_α::Tprog      # tangential boundary condition for the ghost-point copy

    # PHYSICAL CONSTANTS
    g::T                    # gravity
    g_scaled::T             # gravity scaled by scale/scale_η
    cD::T                   # quadratic bottom friction - incl grid spacing
    rD::T                   # linear bottom friction - incl grid spacing
    γ::T                    # frequency of interface relaxation
    cSmag::T                # Smagorinsky constant
    νB::T                   # biharmonic diffusion coefficient
    τSST::T                 # tracer restoring timescale
    jSST::T                 # tracer consumption timescale
    ωFη::Float64            # frequency [1/s] of seasonal surface forcing incl 2π
    ωFx::Float64            # frequency [1/s] of seasonal wind x incl 2π
    ωFy::Float64            # frequency [1/2] of seasonal wind y incl 2π

    # SCALING
    scale::T                # multiplicative constant for low-precision arithmetics used for u,v
    scale_η::T              # multiplicative constant for η
    scale_inv::T            # and its inverse
    scale_sst::T            # scale for sst
end

"""Generator function for the mutable struct Constants."""
function Constants{T,Tprog}(P::Parameter,G::Grid) where {T<:AbstractFloat,Tprog<:AbstractFloat}

    # Runge-Kutta 2nd/3rd/4th order coefficients including time step Δt and grid spacing Δ
    # a are the coefficents to sum the rhs on the fly, such that sum=1
    # b are the coefficents for the update that used for a new evaluation of the RHS
    if P.RKo == 2     # Heun's method
        RKaΔt = Tprog.([1/2,1/2]*G.dtint/G.Δ)
        RKbΔt = Tprog.([1]*G.dtint/G.Δ)
    elseif P.RKo == 3     # version 2 / Heun's 3rd order
        RKaΔt = Tprog.([1/4,0.,3/4]*G.dtint/G.Δ)
        RKbΔt = Tprog.([1/3,2/3]*G.dtint/G.Δ)
    elseif P.RKo == 4
        RKaΔt = Tprog.([1/6,1/3,1/3,1/6]*G.dtint/G.Δ)
        RKbΔt = Tprog.([.5,.5,1.]*G.dtint/G.Δ)
    end

    # Δt/(s-1) for SSPRK2
    Δt_Δs = convert(Tprog,G.dtint/G.Δ/(P.RKs-1))

    # time step and half the time step including the grid spacing as this is not included in the RHS
    Δt_Δ = convert(Tprog,G.dtint/G.Δ)
    Δt_Δ_half = convert(Tprog,G.dtint/G.Δ/2)

    # coefficients for SSPRK3
    SSPRK3c = SSPRK3coeff{Tprog}(P,Δt_Δ)

    # BOUNDARY CONDITIONS AND PHYSICS
    one_minus_α = convert(Tprog,1-P.α)      # for the ghost point copy/tangential boundary conditions
    g = convert(T,P.g)                      # gravity - for Bernoulli potential
    g_scaled = convert(T,P.g*P.scale/P.scale_η)

    # BOTTOM FRICTION COEFFICIENTS
    # incl grid spacing Δ for non-dimensional gradients
    # include scale for quadratic cD only to unscale the scale^2 in u^2
    cD = convert(T,-G.Δ*P.cD/P.scale)     # quadratic drag [m]
    rD = convert(T,-G.Δ/(P.τD*24*3600))   # linear drag [m/s]

    # INTERFACE RELAXATION FREQUENCY
    # incl grid spacing Δ for non-dimensional gradients
    γ = convert(T,G.Δ/(P.t_relax*3600*24))    # [m/s]

    # BIHARMONIC DIFFUSION
    # undo scaling here as smagorinksy diffusion contains scale^2 due to ~u^2
    cSmag = convert(T,-P.cSmag/P.scale)             # Smagorinsky coefficient

    # νB ~ Δ²νA, biharmonic νB is scaled with Δ² from harmonic νA
    # νA ~ (Δ²/(30km)^2)*νA0, i.e. νA is scaled with Δ² relative to νA0=500 at 30km.
    # dimensionless operators: scaled νB* = Δ⁻³νB
    νB = convert(T,-P.νA0*(G.Δ/30000^2))             # Δ² scaling based on 500m^s/s at Δ=30km

    # TRACER ADVECTION
    τSST = convert(T,G.dtadvint/(P.τSST*3600*24))   # tracer restoring [1]
    jSST = convert(T,G.dtadvint/(P.jSST*3600*24))   # tracer consumption [1]

    @unpack tracer_relaxation, tracer_consumption = P
    τSST = tracer_relaxation ? τSST : zero(T)       # set zero as τ,j will be added   
    jSST = tracer_consumption ? jSST : zero(T)      # and executed in one loop

    # TIME DEPENDENT FORCING
    ωFη = -2π*P.ωFη/24/365.25/3600
    ωFx = 2π*P.ωFx/24/365.25/3600
    ωFy = 2π*P.ωFy/24/365.25/3600

    # SCALING
    scale = convert(T,P.scale)
    scale_η = convert(T,P.scale_η)
    scale_inv = convert(T,1/P.scale)
    scale_sst = convert(T,P.scale_sst)

    return Constants{T,Tprog}(  RKaΔt,RKbΔt,Δt_Δs,Δt_Δ,Δt_Δ_half,
                                SSPRK3c,one_minus_α,
                                g,g_scaled,cD,rD,γ,cSmag,νB,τSST,jSST,
                                ωFη,ωFx,ωFy,
                                scale,scale_η,scale_inv,scale_sst)
end