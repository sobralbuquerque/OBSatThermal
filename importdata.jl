using DelimitedFiles, Interpolations

moving_average(vs, n) = [sum(@view vs[i:(i+n-1)]) / n for i = 1:(length(vs)-(n-1))]

function complete_vector!(avgd, vec)
    n_pad = (length(vec) - length(avgd)) ÷ 2

    prepend!(avgd, repeat([vec[begin]], n_pad))
    append!(avgd, repeat([vec[end]], n_pad))
end

complete_moving_average(vs, n) = complete_vector!(moving_average(vs, n), vs)


function read_flight_values(filename::AbstractString)
    data = readdlm(filename, ',', Float64, header = true)[1]
    t, v, T, h = eachcol(data)
    t .-= t[begin]
    t *= u"s"
    v *= u"m/s"
    T = T * u"°C" .|> u"K"
    h *= u"m"
    v = complete_moving_average(v, 35)
    t0, tf = t[begin], t[end]
    dt = (tf - t0) / (length(t) - 1)
    nodes = (t0:dt:tf,)
    return t, v, T, h, nodes
end

data_filename = "flight.csv"
t, v, T, h, t_nodes = read_flight_values(data_filename)
const flight_time = t
const flight_speed = interpolate(t_nodes, v, Gridded(Linear()))
const flight_temperature = interpolate(t_nodes, T, Gridded(Linear()))
const flight_altitude = interpolate(t_nodes, h, Gridded(Linear()))
