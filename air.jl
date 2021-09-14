import CoolProp

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

function pressure(h, T)
    P₀ = 1u"atm"
    M = 0.0289644u"kg/mol"
    R = Unitful.R
    P₀ * exp(-g * M * h / (R * T))
end



function setvalues!(air::Air, t)
    air.T = flight_temperature(t)
    air.u = flight_speed(t)
    h = flight_altitude(t)
    air.P = pressure(h, air.T)
end