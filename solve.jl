include("thermalsystem.jl")
using UnitfulRecipes, Plots, DifferentialEquations
plotly()


T0 = 25.0u"°C" |> u"K"
k_iso = 0.032u"W/m/K"
L_iso = 50.0u"mm"
R_wall = L_iso / k_iso
Cp_battery = 800 * 45.8e-3u"J/K"
L_battery = [9.2, 21.5, 27.3]u"mm"
Q_battery = 50.0u"mW"
T_ext = -50.0u"°C" |> u"K"



function f(T, p, t)
    L_battery, Q_battery = p
    thermal_system(T, T_ext, P_air, R_wall, Cp_battery, L_battery, Q_battery)
end
p = (L_battery, Q_battery)

t = (0.0:110.0)u"minute"
tspan = (t[begin], t[end])

prob = ODEProblem(f, T0, tspan, p)
sol = solve(prob, Tsit5(), saveat = t)
plot(sol.t, sol.u, yunit = u"°C", xlabel="Time", ylabel="Temperature")
