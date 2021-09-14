using UnitfulRecipes, Plots, DifferentialEquations
include("thermalsystem.jl")
plotly()
theme(:mute)

T0 = 25.0u"°C" |> u"K"
k_iso = 0.032u"W/m/K"
Cp_battery = 800 * 45.8e-3u"J/K"
L_battery = [9.2, 21.5, 27.3]u"mm"
Q_battery = 2u"W"

t = flight_time .|> u"minute"
tspan = (t[begin], t[end])

air = Air(T0)
isopor = Wall(T0, 1u"mm" / k_iso)
interface = Interface(air, isopor)

function f(u, p, t)
    T = u
    thermal_system!(interface, T, t, Cp_battery, L_battery, Q_battery)
end

function reset!(i::Interface, R)
    i.air.T = T0
    i.air.u = 0u"m/s"
    i.air.P = 1u"atm"
    i.wall.T = T0
    i.wall.R = R
end

begin
    plot()
    for L ∈ (0.0:2.5:5.0)u"mm"
        reset!(interface, L / k_iso)
        prob = ODEProblem(f, T0, tspan)
        sol = solve(prob, Tsit5(), saveat = t)
        plot!(sol.t, sol.u, label = "L = $L", yunit = u"°C")
    end

    plot!(
        t,
        flight_temperature.(t),
        label = "Air",
        linestyle = :dash,
        color = :black,
        xlabel = "Time",
        ylabel = "Temperature",
    )
end
