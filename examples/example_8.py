import ctypes
import numpy as np
from agdeblend import deblend

n_traces = 128
n_times = [512, 768]
volumes = [0, 1]
n_dims = [2] * 2
window_shapes = [[64, 256]] * 2
coords = [[0]] * 2
shapes = [[n_traces, n_times_v] for n_times_v in n_times]
shottimes = [
    np.fromfile("data_3_shottimes_{}.bin".format(v), ctypes.c_long)
    for v in range(2)
]
channels = [
    np.fromfile("data_3_channels_{}.bin".format(v), ctypes.c_int)
    for v in range(2)
]
trace_types = [
    np.fromfile("data_3_trace_types_{}.bin".format(v), ctypes.c_int)
    for v in range(2)
]
data = [
    np.fromfile(
        "data_3_blended_data_{}.bin".format(v), ctypes.c_float
    ).reshape(shapes[v])
    for v in range(2)
]
initial_factor = 1.0
n_its = 2500
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

for volume_idx in range(2):
    filename = "out/data_3_deblended_data_{}.bin".format(volume_idx)
    data[volume_idx].tofile(filename)
