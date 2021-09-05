import ctypes
import numpy as np
from agdeblend import deblend

n_traces = 128
n_times = 512
volumes = [0]
n_dims = [2]
window_shapes = [[32, 256]]
coords = [[0]]
shapes = [[n_traces, n_times]]
shottimes = [np.fromfile("data_1_shottimes.bin", ctypes.c_long)]
channels = [np.fromfile("data_1_channels.bin", ctypes.c_int)]
trace_types = [np.fromfile("data_1_trace_types.bin", ctypes.c_int)]
data = [
    np.fromfile("data_1_blended_data.bin", ctypes.c_float).reshape(shapes[0])
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

data[0].tofile("out/data_1_deblended_data.bin")
