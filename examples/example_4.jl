include("../wrappers/julia/agdeblend.jl")
using .AGDeblend

n_traces = 128
n_times = 512
volumes = [0]
n_dims = [2]
window_shapes = [[256, 32]]
coords = [[0]]
shapes = [[n_times, n_traces]]
shottimes = [zeros(Clong, n_traces)]
read!("data_1_shottimes.bin", shottimes[1])
channels = [zeros(Cint, n_traces)]
read!("data_1_channels.bin", channels[1])
trace_types = [zeros(Cint, n_traces)]
read!("data_1_trace_types.bin", trace_types[1])
data = [zeros(Cfloat, (n_times, n_traces))]
read!("data_1_blended_data.bin", data[1])
initial_factor = 1.0
n_its = 1000
print_freq = -1

deblend(volumes, window_shapes, coords, shottimes, channels, trace_types,
        initial_factor, n_its, print_freq, data)

write("out/data_1_deblended_data.bin", data[1])
