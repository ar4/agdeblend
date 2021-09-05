include("../wrappers/julia/agdeblend.jl")
using .AGDeblend

n_traces = 128
n_times = [512, 768]
volumes = [0, 1]
n_dims = [2, 2]
window_shapes = [[256, 64], [256, 64]]
coords = [[0], [0]]
shottimes = [zeros(Clong, n_traces), zeros(Clong, n_traces)]
read!("data_3_shottimes_0.bin", shottimes[1])
read!("data_3_shottimes_1.bin", shottimes[2])
channels = [zeros(Cint, n_traces), zeros(Cint, n_traces)]
read!("data_3_channels_0.bin", channels[1])
read!("data_3_channels_1.bin", channels[2])
trace_types = [zeros(Cint, n_traces), zeros(Cint, n_traces)]
read!("data_3_trace_types_0.bin", trace_types[1])
read!("data_3_trace_types_1.bin", trace_types[2])
data = [zeros(Cfloat, (n_times[1], n_traces)),
        zeros(Cfloat, (n_times[2], n_traces))]
read!("data_3_blended_data_0.bin", data[1])
read!("data_3_blended_data_1.bin", data[2])
initial_factor = 1.0
n_its = 2500
print_freq = -1

deblend(volumes, window_shapes, coords, shottimes, channels, trace_types,
        initial_factor, n_its, print_freq, data)

write("out/data_3_deblended_data_0.bin", data[1])
write("out/data_3_deblended_data_1.bin", data[2])
