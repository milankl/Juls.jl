"""Transit function to call either the rhs_linear or the rhs_nonlinear."""
function rhs!(  u::Array{T,2},
                v::Array{T,2},
                η::Array{T,2},
                Diag::DiagnosticVars{T,Tprog},
                S::ModelSetup{T,Tprog},
                t::Int) where {T,Tprog}

    @unpack dynamics = S.parameters

    if dynamics == "linear"
        rhs_linear!(u,v,η,Diag,S,t)
    else
        rhs_nonlinear!(u,v,η,Diag,S,t)
    end
end

"""Tendencies du,dv of

        ∂u/∂t = qhv - ∂(1/2*(u²+v²) + gη)/∂x + Fx
        ∂v/∂t = -qhu - ∂(1/2*(u²+v²) + gη)/∂y + Fy

the nonlinear shallow water equations."""
function rhs_nonlinear!(u::AbstractMatrix,
                        v::AbstractMatrix,
                        η::AbstractMatrix,
                        Diag::DiagnosticVars,
                        S::ModelSetup,
                        t::Int)

    @unpack h,h_u,h_v,U,V = Diag.VolumeFluxes
    @unpack H = S.forcing
    @unpack ep = S.grid

    UVfluxes!(u,v,η,Diag,S)              # U,V needed for PV advection and in the continuity equation
    if S.grid.nstep_advcor == 0              # evaluate every RK substep
        advection_coriolis!(u,v,η,Diag,S)    # PV and non-linear Bernoulli terms
    end
    PVadvection!(Diag,S)                 # advect the PV with U,V
    
    # Bernoulli potential - recalculate for new η, KEu,KEv are only updated in advection_coriolis
    @unpack p,KEu,KEv,dpdx,dpdy = Diag.Bernoulli
    @unpack g_scaled,scale_inv = S.constants
    bernoulli!(p,KEu,KEv,η,g_scaled,ep,scale_inv)
    ∂x!(dpdx,p)
    ∂y!(dpdy,p)

    # adding the terms
    momentum_u!(Diag,S,t)
    momentum_v!(Diag,S,t)
end

"""Tendencies du,dv of

        ∂u/∂t = fv - g∂η/∂x + Fx
        ∂v/∂t = -fu - g∂η/∂y + Fy

the linear shallow water equations."""
function rhs_linear!(   u::AbstractMatrix,
                        v::AbstractMatrix,
                        η::AbstractMatrix,
                        Diag::DiagnosticVars,
                        S::ModelSetup,
                        t::Int)

    @unpack h,h_u,h_v,U,V,dUdx,dVdy = Diag.VolumeFluxes
    @unpack g,scale,scale_η = S.constants
    @unpack ep = S.grid

    # Pressure gradient
    @unpack dpdx,dpdy = Diag.Bernoulli
    g_scale = g*scale/scale_η
    ∂x!(dpdx,g_scale*η)
    ∂y!(dpdy,g_scale*η)

    # Coriolis force
    @unpack qhv,qhu,v_u,u_v = Diag.Vorticity
    @unpack f_u,f_v = S.grid
    Ixy!(v_u,v)
    Ixy!(u_v,u)
    fv!(qhv,f_u,v_u,ep)
    fu!(qhu,f_v,u_v,ep)

    # adding the terms
    momentum_u!(Diag,S,t)
    momentum_v!(Diag,S,t)
end

""" Update advective and Coriolis tendencies."""
function advection_coriolis!(   u::Array{T,2},
                                v::Array{T,2},
                                η::Array{T,2},
                                Diag::DiagnosticVars{T,Tprog},
                                S::ModelSetup{T,Tprog}) where {T,Tprog}

    @unpack h = Diag.VolumeFluxes
    @unpack h_q,dvdx,dudy = Diag.Vorticity
    @unpack u²,v²,KEu,KEv = Diag.Bernoulli

    Ixy!(h_q,h)

    # off-diagonals of stress tensor ∇(u,v)
    ∂x!(dvdx,v)
    ∂y!(dudy,u)

    # non-linear part of the Bernoulli potential
    speed!(u²,v²,u,v)
    Ix!(KEu,u²)
    Iy!(KEv,v²)

    # Potential vorticity update
    PV!(Diag,S)

    @unpack q = Diag.Vorticity          # now load what just has been calculated

    # Linear combinations of the potential vorticity q
    if S.parameters.adv_scheme == "Sadourny"
        @unpack q_u,q_v = Diag.Vorticity
        Iy!(q_u,q)
        Ix!(q_v,q)
    elseif S.parameters.adv_scheme == "ArakawaHsu"
        @unpack qα,qβ,qγ,qδ = Diag.ArakawaHsu
        AHα!(qα,q)
        AHβ!(qβ,q)
        AHγ!(qγ,q)
        AHδ!(qδ,q)
    end
end

"""Layer thickness h obtained by adding sea surface height η to bottom height H."""
function thickness!(h::AbstractMatrix,η::AbstractMatrix,H::AbstractMatrix)
    m,n = size(h)
    @boundscheck (m,n) == size(η) || throw(BoundsError())
    @boundscheck (m,n) == size(H) || throw(BoundsError())

    @inbounds for i in eachindex(η)
        h[i] = η[i] + H[i]
    end
end

"""Squared velocities u²,v²."""
function speed!(u²::AbstractMatrix,
                v²::AbstractMatrix,
                u::AbstractMatrix,
                v::AbstractMatrix)
    m,n = size(u²)
    @boundscheck (m,n) == size(u) || throw(BoundsError())

    @inbounds for i in eachindex(u)
        u²[i] = u[i]^2
    end

    m,n = size(v²)
    @boundscheck (m,n) == size(v) || throw(BoundsError())

    @inbounds for i in eachindex(v)
        v²[i] = v[i]^2
    end
end

"""Bernoulli potential p = 1/2*(u² + v²) + gη."""
function bernoulli!(p::Array{T,2},
                    KEu::Array{T,2},
                    KEv::Array{T,2},
                    η::Array{T,2},
                    g::T,
                    ep::Int,
                    scale_inv::T) where {T<:AbstractFloat}

    m,n = size(p)
    @boundscheck (m+ep,n+2) == size(KEu) || throw(BoundsError())
    @boundscheck (m+2,n) == size(KEv) || throw(BoundsError())
    @boundscheck (m,n) == size(η) || throw(BoundsError())

    one_half_scale_inv = convert(T,0.5)*scale_inv

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            p[i,j] = one_half_scale_inv*(KEu[i+ep,j+1] + KEv[i+1,j]) + g*η[i,j]
        end
    end
end

"""Coriolis term f*v. """
function fv!(   qhv::AbstractMatrix,
                f_u::AbstractMatrix,
                v_u::AbstractMatrix,
                ep::Int)

    m,n = size(qhv)
    @boundscheck (m,n) == size(f_u) || throw(BoundsError())
    @boundscheck (m+4-ep,n+2) == size(v_u) || throw(BoundsError())
    @boundscheck ep == 1 || ep == 0 || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            qhv[i,j] = f_u[i,j]*v_u[i+2-ep,j+1]
        end
    end
end

"""Coriolis term f*u. """
function fu!(   qhu::AbstractMatrix,
                f_v::AbstractMatrix,
                u_v::AbstractMatrix,
                ep::Int)

    m,n = size(qhu)
    @boundscheck (m,n) == size(f_v) || throw(BoundsError())
    @boundscheck (m+2+ep,n+4) == size(u_v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            qhu[i,j] = f_v[i,j]*u_v[i+1+ep,j+2]
        end
    end
end

"""Sum up the tendencies of the non-diffusive right-hand side for the u-component."""
function momentum_u!(   Diag::DiagnosticVars{T,Tprog},
                        S::ModelSetup,
                        t::Int) where {T,Tprog}

    @unpack du = Diag.Tendencies
    @unpack qhv = Diag.Vorticity
    @unpack dpdx = Diag.Bernoulli
    @unpack Fx = S.forcing
    @unpack ep,halo = S.grid

    m,n = size(du) .- (2halo,2halo)     # cut off the halo
    @boundscheck (m,n) == size(qhv) || throw(BoundsError())
    @boundscheck (m+2-ep,n+2) == size(dpdx) || throw(BoundsError())
    @boundscheck (m,n) == size(Fx) || throw(BoundsError())

    if S.parameters.seasonal_wind_x
        @unpack ωFx = S.constants
        Fxt = Ftime(T,t,ωFx)
    else
        Fxt = one(T)
    end

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            du[i+2,j+2] = (Tprog(qhv[i,j]) - Tprog(dpdx[i+1-ep,j+1])) + Tprog(Fxt*Fx[i,j])
        end
    end
end

"""Sum up the tendencies of the non-diffusive right-hand side for the v-component."""
function momentum_v!(   Diag::DiagnosticVars{T,Tprog},
                        S::ModelSetup,
                        t::Int) where {T,Tprog}

    @unpack dv = Diag.Tendencies
    @unpack qhu = Diag.Vorticity
    @unpack dpdy = Diag.Bernoulli
    @unpack Fy = S.forcing
    @unpack halo = S.grid

    m,n = size(dv) .- (2halo,2halo)     # cut off the halo
    @boundscheck (m,n) == size(qhu) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(dpdy) || throw(BoundsError())

    if S.parameters.seasonal_wind_y
        @unpack ωFy = S.constants
        Fyt = Ftime(T,t,ωFy)
    else
        Fyt = one(T)
    end

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
             dv[i+2,j+2] = -(Tprog(qhu[i,j]) + Tprog(dpdy[i+1,j+1])) + Tprog(Fyt*Fy[i,j])
        end
    end
end

"""Zonal mass flux U = uh."""
function Uflux!(U::AbstractMatrix{T},
                u::AbstractMatrix{T},
                h_u::AbstractMatrix{T},
                ep::Int,
                scale_inv::T) where {T<:AbstractFloat}

    m,n = size(U)
    @boundscheck (m,n) == size(h_u) || throw(BoundsError())
    @boundscheck (m+2+ep,n+2) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            U[i,j] = u[1+ep+i,1+j]*h_u[i,j]*scale_inv
        end
    end
end

"""Meridional mass flux V = vh."""
function Vflux!(V::AbstractMatrix{T},
                v::AbstractMatrix{T},
                h_v::AbstractMatrix{T},
                scale_inv::T) where {T<:AbstractFloat}

    m,n = size(V)
    @boundscheck (m,n) == size(h_v) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            V[i,j] = v[i+1,j+1]*h_v[i,j]*scale_inv
        end
    end
end

"""Calculate the mass/volume fluxes U,V from u,v,η."""
function UVfluxes!( u::AbstractMatrix,
                    v::AbstractMatrix,
                    η::AbstractMatrix,
                    Diag::DiagnosticVars,
                    S::ModelSetup)
    
    @unpack h,h_u,h_v,U,V = Diag.VolumeFluxes
    @unpack H = S.forcing
    @unpack ep = S.grid
    @unpack scale_inv = S.constants

    thickness!(h,η,H)
    Ix!(h_u,h)
    Iy!(h_v,h)

    # mass or volume flux U,V = uh,vh
    Uflux!(U,u,h_u,ep,scale_inv)
    Vflux!(V,v,h_v,scale_inv)
end