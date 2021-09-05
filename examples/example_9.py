import ctypes
import numpy as np
from agdeblend import deblend, AGDMissing

n_times = 512
volumes = []
n_dims = [3]
window_shapes = [[16, 16, 256]]
coords = []
space_shape_all = [32, 40]  # 40 = 2 * 24 - 8
file_0_range = (slice(None, 16), slice(None, 24))
file_1_range = (slice(16, 32), slice(16, 40))
file_shape = [16, 24]
patch_ranges_x = [slice(0, 16), slice(8, 24), slice(16, 32)]
patch_ranges_y = [slice(0, 16), slice(8, 32), slice(24, 40)]
patch_ranges = []
for ix in range(3):
    for iy in range(3):
        if (ix == 2 and iy == 0) or (ix == 0 and iy == 2):
            continue
        patch_ranges.append((patch_ranges_x[ix], patch_ranges_y[iy]))
        coords.append([ix, iy])
        volumes.append(0)

shottimes_all = np.zeros(space_shape_all, ctypes.c_long)
shottimes_all[file_0_range] = np.fromfile(
    "data_4_shottimes_0.bin", ctypes.c_long
).reshape(file_shape)
shottimes_all[file_1_range] = np.fromfile(
    "data_4_shottimes_1.bin", ctypes.c_long
).reshape(file_shape)
shottimes = [shottimes_all[patch_range].copy() for patch_range in patch_ranges]
del shottimes_all

channels_all = np.zeros(space_shape_all, ctypes.c_int)
channels_all[file_0_range] = np.fromfile(
    "data_4_channels_0.bin", ctypes.c_int
).reshape(file_shape)
channels_all[file_1_range] = np.fromfile(
    "data_4_channels_1.bin", ctypes.c_int
).reshape(file_shape)
channels = [channels_all[patch_range].copy() for patch_range in patch_ranges]
del channels_all

trace_types_all = AGDMissing * np.ones(space_shape_all, ctypes.c_int)
trace_types_all[file_0_range] = np.fromfile(
    "data_4_trace_types_0.bin", ctypes.c_int
).reshape(file_shape)
trace_types_all[file_1_range] = np.fromfile(
    "data_4_trace_types_1.bin", ctypes.c_int
).reshape(file_shape)
trace_types = [
    trace_types_all[patch_range].copy() for patch_range in patch_ranges
]
del trace_types_all

data_all = np.zeros(space_shape_all + [n_times], ctypes.c_float)
data_all[file_0_range] = np.fromfile(
    "data_4_blended_data_0.bin", ctypes.c_float
).reshape(file_shape + [n_times])
data_all[file_1_range] = np.fromfile(
    "data_4_blended_data_1.bin", ctypes.c_float
).reshape(file_shape + [n_times])
data = [data_all[patch_range].copy() for patch_range in patch_ranges]
del data_all

initial_factor = 1.0
n_its = 1000
print_freq = -1

data = deblend(
    volumes,
    window_shapes,
    coords,
    shottimes,
    channels,
    trace_types,
    initial_factor,
    n_its,
    print_freq,
    data,
)

data_all = np.zeros(space_shape_all + [n_times], ctypes.c_float)

for patch_range, patch_data in zip(patch_ranges, data):
    data_all[patch_range] += patch_data

data_all[file_0_range].tofile("out/data_4_deblended_data_0.bin")
data_all[file_1_range].tofile("out/data_4_deblended_data_1.bin")
