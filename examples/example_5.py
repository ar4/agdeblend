import ctypes
import numpy as np
from agdeblend import deblend

n_traces = [80, 64]
n_times = 512
volumes = [0, 0]
n_dims = [2]
window_shapes = [[32, 256]]
coords = [[0], [1]]
shapes = [[n_traces[0], n_times], [n_traces[1], n_times]]
shottimes = [
    np.fromfile("data_1_shottimes.bin", ctypes.c_long, count=n_traces[0]),
    np.fromfile(
        "data_1_shottimes.bin",
        ctypes.c_long,
        offset=64 * ctypes.sizeof(ctypes.c_long),
        count=n_traces[1],
    ),
]
channels = [
    np.fromfile("data_1_channels.bin", ctypes.c_int, count=n_traces[0]),
    np.fromfile(
        "data_1_channels.bin",
        ctypes.c_int,
        offset=64 * ctypes.sizeof(ctypes.c_int),
        count=n_traces[1],
    ),
]
trace_types = [
    np.fromfile("data_1_trace_types.bin", ctypes.c_int, count=n_traces[0]),
    np.fromfile(
        "data_1_trace_types.bin",
        ctypes.c_int,
        offset=64 * ctypes.sizeof(ctypes.c_int),
        count=n_traces[1],
    ),
]
data = [
    np.fromfile(
        "data_1_blended_data.bin", ctypes.c_float, count=n_traces[0] * n_times
    ).reshape(shapes[0]),
    np.fromfile(
        "data_1_blended_data.bin",
        ctypes.c_float,
        offset=64 * n_times * ctypes.sizeof(ctypes.c_float),
        count=n_traces[1] * n_times,
    ).reshape(shapes[1]),
]
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

n_overlap = window_shapes[0][0] // 2
data[0][n_traces[0] - n_overlap :] += data[1][:n_overlap]

with open("out/data_1_deblended_data.bin", mode="w") as fid:
    data[0].tofile(fid)
    data[1][n_overlap:].tofile(fid)
