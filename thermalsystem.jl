using StaticArrays, Unitful

const g = 9.81u"m/s^2"


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


function parameters(air::Air, T_wall, L)
    dT = abs(air.T - T_wall)
    ρ = density(air)
    μ = viscosity(air)
    Pr = prandtl(air)
    k = conductivity(air)
    v = μ / ρ # dynamic viscosity
    α = v / Pr # thermal diffusivity
    β = 1 / air.T
    Ra = (g * β * dT * L^3) / (α * v)
    (k, Pr, NoUnits(Ra))
end

function setvalues!(air::Air, data, t)
    air.T = data.temperature(t)
    air.u = data.velocity(t)
    h = data.altitude(t)
    air.P = pressure(h, air.T)
end

function forced_convection(air::Air, T_wall, L)
    k, Pr, Ra = parameters(air, T_wall, L)

    Re_laminar = 4000
    Re_transient = 8.7e5

    Re_L = (L * air.u) / sym

    if Re_L < Re_transient && Pre > 0.6
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


function convective_vertical(air::Air, T_wall, L)
    k, Pr, Ra = parameters(air, T_wall, L)

    Nu = 0.68 + 0.67 * (Ra^0.25) * (1 + (0.492 / Pr)^(9 / 16)^(-4 / 9))

    (k * Nu) / L
end

function convective_horizontal(air::T, T_wall, L, A, P)
    k, Pr, Ra = parameters(air, T_wall, L)

    if (0.24 < Pr < 2000) && Ra < 1e7
        Nu = 0.14 * Ra^(1 / 3) * (1 + 0.0107Pr) / (1 + 0.01Pr)
    else
        Ra = rayleigh_air(T_med, dT, A / P)
        Nu = (0.56Ra^0.25) / (1 + (0.492 / Pr)^(9 / 16))^(4 / 9)
    end

    (k * Nu) / L
end

R_equiv(R, h_rad, h) = R + 1 / h_rad + 1 / h

function UA_battery(air::Air, T_wall, R_wall, h_rad, L1, L2, L3)
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
    P2 = 2(L2 + L3)

    h_vertical = convective_vertical(air, T_wall, L1)
    h_horizontal = convective_horizontal(air, T_wall, L2, A2, P2)

    R_h = R_equiv(R_wall, h_rad, h_horizontal)
    R_v1 = R_equiv(R_wall, h_rad, h_vertical)
    R_v2 = R_equiv(R_wall, h_rad, h_vertical)
    G_eq = (2 / R_h + 2 / R_v1 + 2 / R_v2)

    At * G_eq
end

function wall_temperatures(air::Air, T_wall, T_battery, R_wall, h_rad, L1, L2, L3)
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
    P2 = 2(L2 + L3)

    UA = UA_battery(air, T_wall, R_wall, h_rad, L1, L2, L3)
    q = (T_battery - T_air) * UA

    h_vertical = convective_vertical(air, T_wall, L1)
    h_horizontal = convective_horizontal(air, T_wall, L2, A2, P2)
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

convective_vertical(air, T0, L_iso)


using Symbolics

begin
    @syms L1 L2 L3
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
end
