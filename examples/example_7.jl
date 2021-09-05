include("../wrappers/julia/agdeblend.jl")
using .AGDeblend

n_traces = 64
n_times = 1024
volumes = [0]
n_dims = [2]
window_shapes = [[256, 32]]
coords = [[0]]
shapes = [[n_times, n_traces]]
shottimes = [zeros(Clong, n_traces)]
read!("data_2_shottimes.bin", shottimes[1])
channels = [zeros(Cint, n_traces)]
read!("data_2_channels.bin", channels[1])
trace_types = [zeros(Cint, n_traces)]
read!("data_2_trace_types.bin", trace_types[1])
data = [zeros(Cfloat, (n_times, n_traces))]
read!("data_2_blended_data.bin", data[1])
wavelets = [zeros(Cfloat, (513))]
read!("data_2_wavelet.bin", wavelets[1])
wavelet_idxs = [zeros(Cint, (n_traces))]
read!("data_2_wavelet_idxs.bin", wavelet_idxs[1])
initial_factor = 1.0
n_its = 1250
print_freq = -1

deblend(volumes, window_shapes, coords, shottimes, channels, trace_types,
        initial_factor, n_its, print_freq, data, wavelet_idxs=wavelet_idxs,
        wavelets=wavelets)

write("out/data_2_deblended_data.bin", data[1])
