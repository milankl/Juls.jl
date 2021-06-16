""" Extends the matrices u,v,η,sst with a halo of ghost points for boundary conditions."""
function add_halo(  u::Array{T,2},
                    v::Array{T,2},
                    η::Array{T,2},
                    sst::Array{T,2},
                    S::ModelSetup) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = S.grid
    @unpack halo,haloη,halosstx,halossty = S.grid

    # Add zeros to satisfy kinematic boundary conditions
    u = cat(zeros(T,halo,nuy),u,zeros(T,halo,nuy),dims=1)
    u = cat(zeros(T,nux+2*halo,halo),u,zeros(T,nux+2*halo,halo),dims=2)

    v = cat(zeros(T,halo,nvy),v,zeros(T,halo,nvy),dims=1)
    v = cat(zeros(T,nvx+2*halo,halo),v,zeros(T,nvx+2*halo,halo),dims=2)

    η = cat(zeros(T,haloη,ny),η,zeros(T,haloη,ny),dims=1)
    η = cat(zeros(T,nx+2*haloη,haloη),η,zeros(T,nx+2*haloη,haloη),dims=2)

    sst = cat(zeros(T,halosstx,ny),sst,zeros(T,halosstx,ny),dims=1)
    sst = cat(zeros(T,nx+2*halosstx,halossty),sst,zeros(T,nx+2*halosstx,halossty),dims=2)

    # SCALING
    @unpack scale,scale_η,scale_sst = S.constants
    u *= scale
    v *= scale
    η *= scale_η
    sst *= scale_sst

    ghost_points!(u,v,η,S)
    ghost_points_sst!(sst,S)
    return u,v,η,sst
end

"""Cut off the halo from the prognostic variables."""
function remove_halo(   u::Array{T,2},
                        v::Array{T,2},
                        η::Array{T,2},
                        sst::Array{T,2},
                        S::ModelSetup) where {T<:AbstractFloat}

    @unpack halo,haloη,halosstx,halossty = S.grid
    @unpack scale_η,scale_inv,scale_sst = S.constants

    # undo scaling as well
    @views ucut = scale_inv*u[halo+1:end-halo,halo+1:end-halo]
    @views vcut = scale_inv*v[halo+1:end-halo,halo+1:end-halo]
    @views ηcut = η[haloη+1:end-haloη,haloη+1:end-haloη]/scale_η
    @views sstcut = sst[halosstx+1:end-halosstx,halossty+1:end-halossty]/scale_sst

    return ucut,vcut,ηcut,sstcut
end

""" Copy ghost points for u from inside to the halo in the nonperiodic case. """
function ghost_points_u_nonperiodic!(u::AbstractMatrix{T},one_minus_α::T) where T
    n,m = size(u)

    # tangential boundary condition
    @inbounds for i in 1:n
        u[i,1] = one_minus_α*u[i,4]
        u[i,2] = one_minus_α*u[i,3]
        u[i,end-1] = one_minus_α*u[i,end-2]
        u[i,end] = one_minus_α*u[i,end-3]
    end
end

""" Copy ghost points for u from inside to the halo in the periodic case. """
function ghost_points_u_periodic!(::Type{Tcomm},u::AbstractMatrix{T},one_minus_α::T) where {Tcomm,T}
    n,m = size(u)

    # periodic bc - mimic MPI communication with (reduced precision) type Tcomm
    @inbounds for j in 1:m
        u[1,j] = T(Tcomm(u[end-3,j]))
        u[2,j] = T(Tcomm(u[end-2,j]))
        u[end-1,j] = T(Tcomm(u[3,j]))
        u[end,j] = T(Tcomm(u[4,j]))
    end

    # tangential bc
    # ghost_points_u_nonperiodic!(u,one_minus_α)
end

""" Copy ghost points for v from inside to the halo in the nonperiodic case. """
function ghost_points_v_nonperiodic!(v::AbstractMatrix{T},one_minus_α::T) where T
    n,m = size(v)

    # tangential boundary condition
    @inbounds for j in 1:m
        v[1,j] = one_minus_α*v[4,j]
        v[2,j] = one_minus_α*v[3,j]
        v[end-1,j] = one_minus_α*v[end-2,j]
        v[end,j] = one_minus_α*v[end-3,j]
    end
end

""" Copy ghost points for v from inside to the halo in the periodic case. """
function ghost_points_v_periodic!(::Type{Tcomm},v::Array{T,2}) where {Tcomm,T}
    n,m = size(v)

    # mimic MPI communication with (reduced precision) type Tcomm
    @inbounds for j in 1:m
        v[1,j] = T(Tcomm(v[end-3,j]))
        v[2,j] = T(Tcomm(v[end-2,j]))
        v[end-1,j] = T(Tcomm(v[3,j]))
        v[end,j] = T(Tcomm(v[4,j]))
    end
end

""" Copy ghost points for η from inside to the halo in the nonperiodic case. """
function ghost_points_η_nonperiodic!(η::AbstractMatrix)
    n,m = size(η)

    # assume no gradients of η across solid boundaries
    # the 4 corner points are copied twice, but it's faster!
    @inbounds for j in 1:m
        η[1,j] = η[2,j]
        η[end,j] = η[end-1,j]
    end

    @inbounds for i in 1:n
        η[i,1] = η[i,2]
        η[i,end] = η[i,end-1]
    end
end

""" Copy ghost points for η from inside to the halo in the periodic case. """
function ghost_points_η_periodic!(::Type{Tcomm},η::Array{T,2}) where {Tcomm,T}
    n,m = size(η)

    # corner points are copied twice, but it's faster!
    # mimic MPI communication with (reduced precision) type Tcomm
    @inbounds for j in 1:m
        η[1,j] = T(Tcomm(η[end-1,j]))
        η[end,j] = T(Tcomm(η[2,j]))
    end

    @inbounds for i in 1:n
        η[i,1] = η[i,2]
        η[i,end] = η[i,end-1]
    end
end

""" Copy ghost points for η from inside to the halo in the nonperiodic case. """
function ghost_points_sst_nonperiodic!(sst::AbstractMatrix,
                                        halosstx::Int,
                                        halossty::Int)
    n,m = size(sst)

    # assume no gradients of η across solid boundaries
    # the 4 corner points are copied twice, but it's faster!
    @inbounds for j in 1:m
        for i ∈ 1:halosstx
            sst[i,j] = sst[halosstx+1,j]
            sst[end-i+1,j] = sst[end-halosstx,j]
        end
    end

    @inbounds for j ∈ 1:halossty
        for i in 1:n
            sst[i,j] = sst[i,halossty+1]
            sst[i,end-j+1] = sst[i,end-halossty]
        end
    end
end

""" Copy ghost points for η from inside to the halo in the periodic case. """
function ghost_points_sst_periodic!(::Type{Tcomm},
                                    sst::Array{T,2},
                                    halosstx::Int,
                                    halossty::Int) where {Tcomm,T}
    n,m = size(sst)

    # corner points are copied twice, but it's faster!
    @inbounds for j in 1:m
        for i ∈ 1:halosstx
            sst[i,j] = T(Tcomm(sst[end-2*halosstx+i,j]))
            sst[end-halosstx+i,j] = T(Tcomm(sst[halosstx+i,j]))
        end
    end

    @inbounds for j ∈ 1:halossty
        for i in 1:n
            sst[i,j] = sst[i,halossty+1]
            sst[i,end-j+1] = sst[i,end-halossty]
        end
    end
end

"""Decide on boundary condition P.bc which ghost point function to execute."""
function ghost_points!( u::AbstractMatrix,
                        v::AbstractMatrix,
                        η::AbstractMatrix,
                        S::ModelSetup)

    @unpack bc,Tcomm = S.parameters
    @unpack one_minus_α = S.constants

    if bc == "periodic"
        ghost_points_u_periodic!(Tcomm,u,one_minus_α)
        ghost_points_v_periodic!(Tcomm,v)
        ghost_points_η_periodic!(Tcomm,η)
    else
        ghost_points_u_nonperiodic!(u,one_minus_α)
        ghost_points_v_nonperiodic!(v,one_minus_α)
        ghost_points_η_nonperiodic!(η)
    end
end

"""Decide on boundary condition P.bc which ghost point function to execute."""
function ghost_points_uv!( u::AbstractMatrix,
                        v::AbstractMatrix,
                        S::ModelSetup)

    @unpack bc,Tcomm = S.parameters
    @unpack one_minus_α = S.constants

    if bc == "periodic"
        ghost_points_u_periodic!(Tcomm,u,one_minus_α)
        ghost_points_v_periodic!(Tcomm,v)
    else
        ghost_points_u_nonperiodic!(u,one_minus_α)
        ghost_points_v_nonperiodic!(v,one_minus_α)
    end
end

"""Decide on boundary condition P.bc which ghost point function to execute."""
function ghost_points_η!(   η::AbstractMatrix,
                            S::ModelSetup)

    @unpack bc, Tcomm = S.parameters

    if bc == "periodic"
        ghost_points_η_periodic!(Tcomm,η)
    else
        ghost_points_η_nonperiodic!(η)
    end
end

"""Decide on boundary condition P.bc which ghost point function to execute."""
function ghost_points_sst!(sst::AbstractMatrix,S::ModelSetup)

    @unpack bc,Tcomm = S.parameters
    @unpack halosstx,halossty = S.grid

    if bc == "periodic"
        ghost_points_sst_periodic!(Tcomm,sst,halosstx,halossty)
    else
        ghost_points_sst_nonperiodic!(sst,halosstx,halossty)
    end
end