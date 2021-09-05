include("../wrappers/julia/agdeblend.jl")
using .AGDeblend

n_traces = 128
n_times = 512
volumes = []
n_dims = [3]
window_shapes = [[256, 16, 16]]
coords = []
space_shape_all = [40, 32]
file_0_range = [1:24, 1:16]
file_1_range = [17:40, 17:32]
file_shape = [24, 16]
patch_ranges_y = [1:16, 9:32, 25:40]
patch_ranges_x = [1:16, 9:24, 17:32]
patch_ranges = []
for ix in 1:3
  for iy in 1:3
    if (ix == 3 && iy == 1) || (ix == 1 && iy == 3)
      continue
    end
    append!(patch_ranges, [[patch_ranges_y[iy], patch_ranges_x[ix]]])
    append!(coords, [[iy, ix]])
    append!(volumes, 0)
  end
end
coords = Array{Array{Integer}}(coords)
volumes = Array{Integer}(volumes)

shottimes_all = zeros(Clong, space_shape_all...)
read!("data_4_shottimes_0.bin", view(shottimes_all, file_0_range...))
read!("data_4_shottimes_1.bin", view(shottimes_all, file_1_range...))
shottimes = [shottimes_all[patch_range...] for patch_range in patch_ranges]
shottimes_all = nothing

channels_all = zeros(Cint, space_shape_all...)
read!("data_4_channels_0.bin", view(channels_all, file_0_range...))
read!("data_4_channels_1.bin", view(channels_all, file_1_range...))
channels = [channels_all[patch_range...] for patch_range in patch_ranges]
channels_all = nothing

trace_types_all = AGDMissing * ones(Cint, space_shape_all...)
read!("data_4_trace_types_0.bin", view(trace_types_all, file_0_range...))
read!("data_4_trace_types_1.bin", view(trace_types_all, file_1_range...))
trace_types = [trace_types_all[patch_range...] for patch_range in patch_ranges]
trace_types_all = nothing

data_all = zeros(Cfloat, n_times, space_shape_all...)
read!("data_4_blended_data_0.bin", view(data_all, :, file_0_range...))
read!("data_4_blended_data_1.bin", view(data_all, :, file_1_range...))
data = [data_all[:, patch_range...] for patch_range in patch_ranges]
data_all = nothing

initial_factor = 1.0
n_its = 1000
print_freq = -1

deblend(volumes, window_shapes, coords, shottimes, channels, trace_types,
        initial_factor, n_its, print_freq, data)

data_all = zeros(Cfloat, n_times, space_shape_all...)

for (patch_range, patch_data) in zip(patch_ranges, data)
  data_all[:, patch_range...] += patch_data
end

write("out/data_4_deblended_data_0.bin", view(data_all, :, file_0_range...))
write("out/data_4_deblended_data_1.bin", view(data_all, :, file_1_range...))
