struct Forcing{T<:AbstractFloat}
    Fx::Array{T,2}
    Fy::Array{T,2}
    H::Array{T,2}
    η_ref::Array{T,2}
    Fη::Array{T,2}
    #sst_ref::Array{T,2}
    #sst_γ::Array{T,2}
end

function Forcing{T}(P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack wind_forcing_x,wind_forcing_y,topography = P

    if wind_forcing_x == "channel"
        Fx,_ = ChannelWind(T,P,G)
    elseif wind_forcing_x == "shear"
        Fx,_ = ShearWind(T,P,G)
    elseif wind_forcing_x == "double_gyre"
        Fx,_ = DoubleGyreWind(T,P,G)
    elseif wind_forcing_x == "constant"
        Fx,_ = ConstantWind(T,P,G)
    elseif wind_forcing_x == "none"
        Fx,_ = NoWind(T,P,G)
    end

    if wind_forcing_y == "channel"
        _,Fy = ChannelWind(T,P,G)
    elseif wind_forcing_y == "shear"
        _,Fy = ShearWind(T,P,G)
    elseif wind_forcing_y == "double_gyre"
        _,Fy = DoubleGyreWind(T,P,G)
    elseif wind_forcing_y == "constant"
        _,Fy = ConstantWind(T,P,G)
    elseif wind_forcing_x == "none"
        _,Fy = NoWind(T,P,G)
    end

    if topography == "ridge"
        H,_ = Ridge(T,P,G)
    elseif topography == "ridges"
        H = Ridges(T,P,G)
    elseif topography == "seamount"
        H = Seamount(T,P,G)
    elseif topography == "flat"
        H = FlatBottom(T,P,G)
    end

    η_ref = InterfaceRelaxation(T,P,G)
    Fη = KelvinPump(T,P,G)

    return Forcing{T}(Fx,Fy,H,η_ref,Fη)
end

"""Returns the constant forcing matrices Fx,Fy that vary only meriodionally/zonally
as a cosine with strongest forcing in the middle and vanishing forcing at boundaries."""
function ChannelWind(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack Δ,x_u,y_u,Lx,Ly = G
    @unpack Fx0,Fy0,H,ρ,scale = P

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (scale*Δ*Fx0/ρ/H)*cos.(π*(yy_u/Ly .- 1/2)).^2
    Fy = (scale*Δ*Fy0/ρ/H)*cos.(π*(xx_u/Lx .- 1/2)).^2

    return T.(Fx),T.(Fy)
end

"""Returns the constant forcing matrices Fx,Fy that vary only meriodionally/zonally
as a hyperbolic tangent with strongest shear in the middle."""
function ShearWind(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack Δ,x_u,y_u,Lx,Ly = G
    @unpack Fx0,Fy0,H,ρ,scale = P

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (scale*Δ*Fx0/ρ/H)*tanh.(2π*(yy_u/Ly .- 1/2))
    Fy = (scale*Δ*Fy0/ρ/H)*tanh.(2π*(xx_u/Lx .- 1/2))

    return T.(Fx),T.(Fy)
end

"""Returns the constant forcing matrices Fx,Fy that vary only meriodionally/zonally
with a superposition of sin & cos for a double gyre circulation.
See Cooper&Zanna 2015 or Kloewer et al 2018."""
function DoubleGyreWind(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack Δ,x_u,y_u,Lx,Ly = G
    @unpack Fx0,Fy0,H,ρ,scale = P

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (scale*Δ*Fx0/ρ/H)*(cos.(2π*(yy_u/Ly .- 1/2)) + 2*sin.(π*(yy_u/Ly .- 1/2)))
    Fy = (scale*Δ*Fy0/ρ/H)*(cos.(2π*(xx_u/Lx .- 1/2)) + 2*sin.(π*(xx_u/Lx .- 1/2)))
    return T.(Fx),T.(Fy)
end

"""Returns constant in in space forcing matrices Fx,Fy."""
function ConstantWind(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack Δ,nux,nuy,nvx,nvy = G
    @unpack Fx0,Fy0,H,ρ,scale = P

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    Fx = T.(scale*Δ*Fx0/ρ/H)*ones(T,nux,nuy)
    Fy = T.(scale*Δ*Fy0/ρ/H)*ones(T,nvx,nvy)

    return Fx,Fy
end

"""Returns constant in in space forcing matrices Fx,Fy."""
function NoWind(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack nux,nuy,nvx,nvy = G

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    Fx = zeros(T,nux,nuy)
    Fy = zeros(T,nvx,nvy)

    return Fx,Fy
end

"""Returns a reference state for Newtonian cooling/surface relaxation shaped as a
hyperbolic tangent to force the continuity equation."""
function InterfaceRelaxation(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack x_T,y_T,Ly,Lx = G
    @unpack η_refh,η_refw = P

    # width of a tangent is defined as the distance between 65% of its minimum value and 65% of the max
    xx_T,yy_T = meshgrid(x_T,y_T)
    η_ref = -(η_refh/2)*tanh.(2π*(Ly/(4*η_refw))*(yy_T/Ly .- 1/2))
    return T.(η_ref)
end

"""Returns a matrix of water depth for the whole domain that contains a
Gaussian seamount in the middle. Water depth, heigth and width of the
seamount are adjusted with the constants H, topofeat_height and topofeat_width."""
function Seamount(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack x_T_halo,y_T_halo,Lx,Ly = G
    @unpack topo_width,topo_height,H = P

    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)
    bumpx = exp.(-((xx_T .- Lx/2).^2)/(2*topo_width^2))
    bumpy = exp.(-((yy_T .- Ly/2).^2)/(2*topo_width^2))

    return T.(P.scale_η*(H .- topo_height*bumpx.*bumpy))
end

"""Returns a matrix of water depth for the whole domain that contains a
meridional Gaussian ridge in the middle. Water depth, heigth and width of the
ridge are adjusted with the constants water_depth, topofeat_height and topofeat_width."""
function Ridge(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack x_T_halo,y_T_halo,Lx,Ly = G
    @unpack topo_width,topo_height,H = P

    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)
    bumpx = exp.(-((xx_T .- Lx/2).^2)/(2*topo_width^2))
    bumpy = exp.(-((yy_T .- Ly/2).^2)/(2*topo_width^2))

    Hx = P.scale_η*(H .- topo_height*bumpx)
    Hy = P.scale_η*(H .- topo_height*bumpy)
    return T.(Hx),T.(Hy)
end

"""Same as Ridge() but for n ridges at various x positions."""
function Ridges(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack x_T_halo,y_T_halo,Lx,Ly = G
    @unpack topo_width,topo_height,H = P
    @unpack topo_ridges_positions = P
    n_ridges = length(topo_ridges_positions)

    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)

    # loop over bumps in x direction
    R = zero(xx_T) .+ H

    for i in 1:n_ridges
        R .-= topo_height*exp.(-((xx_T .- topo_ridges_positions[i]*Lx).^2)/(2*topo_width^2))
    end

    return T.(P.scale_η*R)
end

"""Returns a matrix of constant water depth H."""
function FlatBottom(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}
    @unpack nx,ny,haloη = G
    @unpack H = P
    return fill(T(P.scale_η*H),(nx+2*haloη,ny+2*haloη))
end

"""Returns Kelvin wave pumping forcing of the continuity equation."""
function KelvinPump(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack x_T,y_T,Lx,Ly = G
    @unpack R,ϕ,Δ = G
    @unpack A,ϕk,wk = P

    _,yy_T = meshgrid(x_T,y_T)

    mϕ = 2π*R/360.          # meters per degree latitude
    y0 = Ly/2 - (ϕ-ϕk)*mϕ   # y[m] for central latitude of pumping

    Fη = P.scale_η*A*Δ*exp.(-(yy_T.-y0).^2/(2*wk^2))

    return T.(Fη)
end

"""Time evolution of forcing."""
function Ftime(::Type{T},t::Int,ω::Real) where {T<:AbstractFloat}
    return convert(T,sin(ω*t))
end