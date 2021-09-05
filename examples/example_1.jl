include("../wrappers/julia/agdeblend.jl")
using .AGDeblend

n_traces = [128]
n_times = [512]
shottimes = [zeros(Clong, n_traces[1])]
read!("data_1_shottimes.bin", shottimes[1])
channels = [zeros(Cint, n_traces[1])]
read!("data_1_channels.bin", channels[1])
trace_types = [zeros(Cint, n_traces[1])]
read!("data_1_trace_types.bin", trace_types[1])
data = [zeros(Cfloat, (n_times[1], n_traces[1]))]
read!("data_1_true_data.bin", data[1])
blend_mode = AGDBlendSum

blended_data = blend(shottimes, channels, trace_types, data, blend_mode)

write("out/data_1_blended_data.bin", blended_data[1])
