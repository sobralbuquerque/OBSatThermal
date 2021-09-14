using StaticArrays, Unitful
using Unitful: σ
const ε = 0.9
const g = 9.81u"m/s^2"

include("air.jl")
include("interface.jl")
include("importdata.jl")


function forced_convection(i::Interface, L)
    k, Pr, v = i.k, i.Pr, i.v
    u = abs(i.air.u)

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

    h = (k * Nu) / L

    h
end


function convective_vertical(i::Interface, L)
    k, Pr = i.k, i.Pr
    Ra = rayleigh(i, L)

    Nu = 0.68 + 0.67 * (Ra^0.25) * (1 + (0.492 / Pr)^(9 / 16)^(-4 / 9))

    h = (k * Nu) / L

    h
end

function convective_horizontal(i::Interface, L, A, P)
    k, Pr = i.k, i.Pr
    Ra = rayleigh(i, L)
    if (0.24 < Pr < 2000) && Ra < 1e7
        Nu = 0.14Ra^(1 / 3) * (1 + 0.0107Pr) / (1 + 0.01Pr)
    else
        Ra = rayleigh(i, A / P)
        Nu = (0.56Ra^0.25) / (1 + (0.492 / Pr)^(9 / 16))^(4 / 9)
    end

    h = (k * Nu) / L

    h
end

R_equiv(R, h_rad, h) = R + 1 / h_rad + 1 / h

function UA_battery(i::Interface, h_rad, L1, L2, L3)
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
    P2 = 2(L2 + L3)

    if abs(i.air.u) > zero(i.air.u)
        h_vertical = forced_convection(i, L1)
    else
        h_vertical = convective_vertical(i, L1)
    end

    h_horizontal = convective_horizontal(i, L2, A2, P2)
    R_wall = i.wall.R

    R_h = R_equiv(R_wall, h_rad, h_horizontal)
    R_v1 = R_equiv(R_wall, h_rad, h_vertical)
    R_v2 = R_equiv(R_wall, h_rad, h_vertical)
    G_eq = (2 / R_h + 2 / R_v1 + 2 / R_v2)

    UA = At * G_eq

    UA
end

function wall_temperatures!(i::Interface, T_battery, h_rad, L1, L2, L3)
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
    P2 = 2(L2 + L3)

    UA = UA_battery(i, h_rad, L1, L2, L3)
    q = (T_battery - i.air.T) * UA



    if abs(i.air.u) > zero(i.air.u)
        h_vertical = forced_convection(i, L1)
    else
        h_vertical = convective_vertical(i, L1)
    end

    h_horizontal = convective_horizontal(i, L2, A2, P2)


    R = i.wall.R

    R_h = R + 1 / h_horizontal
    R_v1 = R + 1 / h_vertical
    R_v2 = R + 1 / h_vertical
    G_eq = 2 / R_h + 2 / R_v1 + 2 / R_v2

    UA = At * G_eq

    i.wall.T = -q * (1 / UA) + i.air.T
end

function radiation!(i::Interface, T_battery, L_battery)
    error = Inf * u"K"
    tol = 1e-8 * u"K"
    h_rad = 4ε * σ * (i.wall.T^2 + i.air.T^2) * (i.wall.T + i.air.T)
    while error > tol
        T_previous = i.wall.T
        h_rad = 4ε * σ * (i.wall.T^2 + i.air.T^2) * (i.wall.T + i.air.T)
        wall_temperatures!(i, T_battery, h_rad, L_battery...)
        error = abs(i.wall.T - T_previous)
    end
    return h_rad
end

function thermal_system!(i::Interface, T, t, Cp_battery, L_battery, Q_battery)
    T_battery = T
    i.wall.T = T
    setvalues!(i.air, t)
    h_rad = radiation!(i, T_battery, L_battery)
    UA = UA_battery(i, h_rad, L_battery...)
    Q = UA * (i.air.T - T_battery) + Q_battery

    Q / Cp_battery
end
