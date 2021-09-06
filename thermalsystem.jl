import CoolProp
using Unitful
using Unitful: σ
const ε = 0.8
const P_air = 1u"atm"
const g = 9.81u"m/s^2"

air_prop(property_name, T, P) = CoolProp.PropsSI(property_name, "T", T, "P", P, "Air")
air_density(T, P = P_air) = air_prop("Dmass", T, P)
air_viscosity(T, P = P_air) = air_prop("viscosity", T, P) * u"Pa*s"
air_conductivity(T, P = P_air) = air_prop("conductivity", T, P)
air_prandtl(T, P = P_air) = air_prop("Prandtl", T, P)

function rayleigh_air(T_air, T_wall, L, P = P_air)
    T̄ = (T_air + T_wall) / 2
    dT = abs(T_air - T_wall)
    ρ = air_density(T̄, P)
    μ = air_viscosity(T̄, P)
    Pr = air_prandtl(T̄, P)
    v = μ / ρ # dynamic viscosity
    α = v / Pr # thermal diffusivity
    β = 1 / (T_air)
    Ra = (g * β * dT * L^3) / (α * v)
    NoUnits(Ra)
end

function convective_vertical(T_air, T_wall, L, P_air = P_air)
    T_med = (T_air + T_wall) / 2
    k = air_conductivity(T_med, P_air)
    Pr = air_prandtl(T_med, P_air)
    Ra = rayleigh_air(T_air, T_wall, L, P_air)

    # Squire-Eckert
    Nu = 0.68 + 0.67 * (Ra^0.25) * (1 + (0.492 / Pr)^(9 / 16)^(-4 / 9))

    h = (k * Nu) / L
    h
end

function convective_horizontal(T_air, T_wall, L, A, P, P_air = P_air)
    T_med = (T_air + T_wall) / 2
    dT = abs(T_air - T_wall)
    k = air_conductivity(T_med, P_air)
    Pr = air_prandtl(T_med, P_air)
    Ra = rayleigh_air(T_air, T_wall, L, P_air)

    if (0.24 < Pr < 2000) && Ra < 1e7
        Nu = 0.14 * Ra^(1 / 3) * (1 + 0.0107Pr) / (1 + 0.01Pr)
    else
        Ra = rayleigh_air(T_med, dT, A / P, P_air)
        Nu = (0.56Ra^0.25) / (1 + (0.492 / Pr)^(9 / 16))^(4 / 9)
    end

    h = (k * Nu) / L
    h
end

function R_equiv(R, h_rad, A, h)
    R_rad = 1 / (h_rad * A)
    R_air = 1 / (h * A)
    if R_air / unit(R_air) == Inf
        R_eq = R_rad
    else
        R_eq = (R_air * R_rad) / (R_air + R_rad)
    end

    R + R_eq * A
end


function UA_battery(T_air, T_wall, R_wall, h_rad, L1, L2, L3, P_air = P_air)
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
    P2 = 2(L2 + L3)

    h_vertical = convective_vertical(T_air, T_wall, L1, P_air)
    h_horizontal = convective_horizontal(T_air, T_wall, L2, A2, P2, P_air)

    R_h = R_equiv(R_wall, h_rad, A1, h_horizontal)
    R_v1 = R_equiv(R_wall, h_rad, A2, h_vertical)
    R_v2 = R_equiv(R_wall, h_rad, A3, h_vertical)
    G_eq = (2 / R_h + 2 / R_v1 + 2 / R_v2)

    At * G_eq
end

function wall_temperatures(T_air, T_wall, T_battery, P_air, R_wall, h_rad, L1, L2, L3)
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
    P2 = 2(L2 + L3)

    R_eq = 1 / (UA_battery(T_air, T_wall, R_wall, h_rad, L1, L2, L3, P_air))

    q = (T_battery - T_air) / R_eq

    h_vertical = convective_vertical(T_air, T_wall, L1, P_air)
    h_horizontal = convective_horizontal(T_air, T_wall, L2, A2, P2, P_air)

    R_h = R_equiv(R_wall, h_rad, A1, h_horizontal)
    R_v1 = R_equiv(R_wall, h_rad, A2, h_vertical)
    R_v2 = R_equiv(R_wall, h_rad, A3, h_vertical)

    G_eq = (2 / R_h + 2 / R_v1 + 2 / R_v2)

    UA = At * G_eq
    T_wall = -q * (1 / UA) + T_air
end

function thermal_system(T, T_air, P_air, R_wall, Cp_battery, L_battery, Q_battery)
    T_battery = T
    T_wall = T
    error = Inf * u"K"
    tol = 1e-8u"K"
    L1, L2, L3 = L_battery
    h_rad = 4ε * σ * (T_wall^2 + T_air^2) * (T_wall + T_air)
    while error > tol
        T_previous = T_wall
        h_rad = 4ε * σ * (T_wall^2 + T_air^2) * (T_wall + T_air)
        T_wall =
            wall_temperatures(T_air, T_wall, T_battery, P_air, R_wall, h_rad, L1, L2, L3)
        error = abs(T_wall - T_previous)
    end

    UA = UA_battery(T_air, T_wall, R_wall, h_rad, L_battery..., P_air)
    dT_battery = (UA * (T_air - T_battery) + Q_battery) / Cp_battery

    dT_battery
end
