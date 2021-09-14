include("thermalsystem.jl")


T0 = 25.0u"°C" |> u"K"
k_iso = 0.032u"W/m/K"
Cp_battery = 800 * 45.8e-3u"J/K"
L_battery = [9.2, 21.5, 27.3]u"mm"
Q_battery = 0.0u"mW"

t = flight_time .|> u"minute"
tspan = (t[begin], t[end])

air = Air(T0)
isopor = Wall(T0 - 1u"K", 1u"mm" / k_iso)
interface = Interface(air, isopor)


radiation!(interface, T0 - 10u"K", L_battery)

interface.wall.T
interface.air.T
interface.T

i = interface

h_rad = 4ε * σ * (i.wall.T^2 + i.air.T^2) * (i.wall.T + i.air.T)
begin
    dT = 100.0u"K"
    T_battery = T0 + dT
    L1, L2, L3 = L_battery
    A1 = L2 * L3
    A2 = L1 * L2
    A3 = L1 * L3
    At = 2(A1 + A2 + A3)
    P2 = 2(L2 + L3)

    UA = UA_battery(i, h_rad, L1, L2, L3)

    ua1 = UA

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
    UA / UA_battery(i, h_rad, L1, L2, L3) |> NoUnits
    q * (1 / UA) + i.air.T
end
