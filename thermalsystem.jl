import CoolProp
using StaticArrays, Unitful
using Unitful: σ

const ε = 0.8
const g = 9.81u"m/s^2"

air_prop(property_name, T, P) = CoolProp.PropsSI(property_name, "T", T, "P", P, "Air")
air_density(T, P = P_air) = air_prop("Dmass", T, P)
air_viscosity(T, P = P_air) = air_prop("viscosity", T, P)
air_viscosity(T::Quantity, P = P_air) = air_prop("viscosity", T, P) * u"Pa*s"
air_conductivity(T, P = P_air) = air_prop("conductivity", T, P)
air_prandtl(T, P = P_air) = air_prop("Prandtl", T, P)

mutable struct Air
    T::Any
    P::Any
    u::Any
    function Air(T, P = 1u"atm", u = 0u"m/s")
        new(T, P, u)
    end
end

Base.:+(a1::Air, a2::Air) = Air(a1.T + a2.T, a1.P + a2.P, a1.u + a2.u)

density(air::Air) = air_density(air.T, air.P)
viscosity(air::Air) = air_viscosity(air.T, air.P)
conductivity(air::Air) = air_conductivity(air.T, air.P)
prandtl(air::Air) = air_prandtl(air.T, air.P)
mutable struct Wall
    T::Any
    L::Any
    R::Any
end

struct Interface
    air::Air
    wall::Wall
end

density(i::Interface) = air_density(i.T, i.P)
viscosity(i::Interface) = air_viscosity(i.T, i.P)
conductivity(i::Interface) = air_conductivity(i.T, i.P)
prandtl(i::Interface) = air_prandtl(i.T, i.P)

function Base.getproperty(i::Interface, sym::Symbol)
    if sym === :T
        # Average temperature
        (i.air.T + i.wall.T) / 2
    elseif sym === :ΔT
        # Temperature gradient
        abs(i.air.T - i.wall.T)
    elseif sym === :P
        # Pressure
        i.air.P
    elseif sym === :k
        conductivity(i)
    elseif sym === :ρ
        density(i)
    elseif sym === :μ
        viscosity(i)
    elseif sym === :Pr
        prandtl(i)
    elseif sym === :v
        # Dynamic viscosity
        i.μ / i.ρ
    elseif sym === :α
        # Thermal diffusivity
        i.v / i.Pr
    elseif sym === :β
        1 / i.air.T
    elseif sym === :Ra
        # Rayleigh number
        local Ra = g * i.β * i.ΔT * i.wall.L^3 / (i.α * i.v)
        NoUnits(Ra)
    else
        getfield(i, sym)
    end
end



function pressure(h, T)
    P₀ = 1u"atm"
    M = 0.0289644u"kg/mol"
    R = Unitful.R
    P₀ * exp(-g * M * h / (R * T))
end


struct Material
    k::Any
    ρ::Any
end

struct Cuboid{T}
    L::SVector{3,T}
    t::T
    material::Material
end

const DEPRON = Material(
    0.035u"W/m/K", # source:
    50u"kg/m^3", # source: 
)
const ISO40 = Material(
    0.032u"W/m/K", # source:
    22u"kg/m^3", # source:
)
const LI_PO = Material(
    1u"W/m/K", # source:
    1u"kg/m^3", # source:
)

battery = Cuboid(SVector(1, 2, 3)u"mm", 10u"mm", LI_PO)

resistence(c::Cuboid) = c.t / c.material.k



function setvalues!(air::Air, t)
    air.T = temperature(t)
    air.u = velocity(t)
    h = altitude(t)
    air.P = pressure(h, air.T)
end

function forced_convection(i::Interface, L)
    k, Pr = i.k, i.Pr
    u, v = i.air.u, i.air.v

    Re_laminar = 4000
    Re_transient = 8.7e5

    Re_L = (L * u) / v

    if Re_L < Re_transient && Pr > 0.6
        Nu = 0.664 * √Re_L * ∛Pr
    else
        c = 0.992log10(Re_laminar)
        Nu1 = 0.037Pr^0.6 * (Re_L^0.8 - Re_transient^0.8) # Turbulent Nusselt Number
        Nu2 = 0.664 * √Re_L * ∛Pr # Laminar Nusselt Number
        Nu3 = (0.0296Re_transient^0.8 * Pr^0.6 - 0.332 * √Re_laminar * ∛Pr) / c # Transition Nusselt Number
        Nu = Nu1 + Nu2 + Nu3
    end

    (k * Nu) / L
end


function convective_vertical(i::Interface, L)
    k, Pr, Ra = i.k, i.Pr, i.Ra

    Nu = 0.68 + 0.67 * (Ra^0.25) * (1 + (0.492 / Pr)^(9 / 16)^(-4 / 9))

    (k * Nu) / L
end

function convective_horizontal(i::Interface, L, A, P)
    k, Pr, Ra = i.k, i.Pr, i.Ra

    if (0.24 < Pr < 2000) && Ra < 1e7
        Nu = 0.14Ra^(1 / 3) * (1 + 0.0107Pr) / (1 + 0.01Pr)
    else
        Ra = rayleigh_air(T_med, dT, A / P)
        Nu = (0.56Ra^0.25) / (1 + (0.492 / Pr)^(9 / 16))^(4 / 9)
    end

    (k * Nu) / L
end

R_equiv(R, h_rad, h) = R + 1 / h_rad + 1 / h

function UA_battery(i::Interface, h_rad, L1, L2, L3)
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
    P2 = 2(L2 + L3)

    if abs(i.air.u) > 0
        h_vertical = forced_convection(i, L1)
    else
        h_vertical = convective_vertical(i, L1)
    end

    h_horizontal = convective_horizontal(i, L2, A2, P2)

    R_h = R_equiv(R_wall, h_rad, h_horizontal)
    R_v1 = R_equiv(R_wall, h_rad, h_vertical)
    R_v2 = R_equiv(R_wall, h_rad, h_vertical)
    G_eq = (2 / R_h + 2 / R_v1 + 2 / R_v2)

    At * G_eq
end

function wall_temperatures(i, T_battery, h_rad, L1, L2, L3)
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
    P2 = 2(L2 + L3)

    UA = UA_battery(i, h_rad, L1, L2, L3)
    q = (T_battery - T_air) * UA

    h_vertical = convective_vertical(i, L1)
    h_horizontal = convective_horizontal(i, L2, A2, P2)
    R_h = R_wall + 1 / h_horizontal
    R_v1 = R_wall + 1 / h_vertical
    R_v2 = R_wall + 1 / h_vertical
    G_eq = 2 / R_h + 2 / R_v1 + 2 / R_v2

    UA = At * G_eq

    -q * (1 / UA) + T_air
end
##

T0 = 25.0u"°C" |> u"K"
k_iso = 0.032u"W/m/K"
L_iso = 50.0u"mm"
R_wall = L_iso / k_iso
Cp_battery = 800 * 45.8e-3u"J/K"
L_battery = [9.2, 21.5, 27.3]u"mm"
Q_battery = 50.0u"mW"
T_ext = -50.0u"°C" |> u"K"
air = Air(T_ext)
wall = Wall(T_ext, L_iso, R_wall)
interface = Interface(air, wall)
convective_vertical(interface, L_iso)
