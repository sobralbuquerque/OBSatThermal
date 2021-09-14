mutable struct Wall
    T::Any
    R::Any
end

struct Interface
    air::Air
    wall::Wall
end

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
    else
        getfield(i, sym)
    end
end

density(i::Interface) = air_density(i.T, i.P)
viscosity(i::Interface) = air_viscosity(i.T, i.P)
conductivity(i::Interface) = air_conductivity(i.T, i.P)
prandtl(i::Interface) = air_prandtl(i.T, i.P)
rayleigh(i::Interface, L) = g * i.β * i.ΔT * L^3 / (i.α * i.v) |> NoUnits