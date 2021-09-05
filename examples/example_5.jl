include("../wrappers/julia/agdeblend.jl")
using .AGDeblend

n_traces = [80, 64]
n_times = 512
volumes = [0, 0]
n_dims = [2]
window_shapes = [[256, 32]]
coords = [[0], [1]]
shottimes = [zeros(Clong, n_traces[1]), zeros(Clong, n_traces[2])]
open("data_1_shottimes.bin", "r") do f
  read!(f, shottimes[1])
  seek(f, 64 * sizeof(Clong))
  read!(f, shottimes[2])
end
channels = [zeros(Cint, n_traces[1]), zeros(Cint, n_traces[2])]
open("data_1_channels.bin", "r") do f
  read!(f, channels[1])
  seek(f, 64 * sizeof(Cint))
  read!(f, channels[2])
end
trace_types = [zeros(Cint, n_traces[1]), zeros(Cint, n_traces[2])]
open("data_1_trace_types.bin", "r") do f
  read!(f, trace_types[1])
  seek(f, 64 * sizeof(Cint))
  read!(f, trace_types[2])
end
data = [zeros(Cfloat, (n_times, n_traces[1])),
        zeros(Cfloat, (n_times, n_traces[2]))]
open("data_1_blended_data.bin", "r") do f
  read!(f, data[1])
  seek(f, 64 * n_times * sizeof(Cfloat))
  read!(f, data[2])
end
initial_factor = 1.0
n_its = 1000
print_freq = -1

deblend(volumes, window_shapes, coords, shottimes, channels, trace_types,
        initial_factor, n_its, print_freq, data)

n_overlap = div(window_shapes[1][2], 2)
data[1][:, n_traces[1] - n_overlap + 1:end,] += data[2][:, 1:n_overlap]

open("out/data_1_deblended_data.bin", "w") do f
  write(f, data[1])
  write(f, data[2][:, n_overlap + 1:end])
end
