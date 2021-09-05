include("../wrappers/julia/agdeblend.jl")
using .AGDeblend

n_traces = [128]
n_times = [512]
n_times_out = [768]
shottimes = [zeros(Clong, n_traces[1])]
read!("data_1_shottimes.bin", shottimes[1])
shottimes_out = [zeros(Clong, n_traces[1])]
read!("data_1_shottimes_out.bin", shottimes_out[1])
channels = [zeros(Cint, n_traces[1])]
read!("data_1_channels.bin", channels[1])
trace_types = [zeros(Cint, n_traces[1])]
read!("data_1_trace_types.bin", trace_types[1])
data = [zeros(Cfloat, (n_times[1], n_traces[1]))]
read!("data_1_blended_data.bin", data[1])
blend_mode = AGDBlendMean
taper_length = 16

blended_data = blend(shottimes, channels, trace_types, data, blend_mode,
                     taper_length=taper_length, n_times_out=n_times_out,
                     shottimes_out=shottimes_out)

write("out/data_1_blended_data_2.bin", blended_data[1])
