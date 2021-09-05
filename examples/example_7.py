import ctypes
import numpy as np
from agdeblend import deblend

n_traces = 64
n_times = 1024
volumes = [0]
n_dims = [2]
window_shapes = [[32, 256]]
coords = [[0]]
shapes = [[n_traces, n_times]]
shottimes = [np.fromfile("data_2_shottimes.bin", ctypes.c_long)]
channels = [np.fromfile("data_2_channels.bin", ctypes.c_int)]
trace_types = [np.fromfile("data_2_trace_types.bin", ctypes.c_int)]
data = [
    np.fromfile("data_2_blended_data.bin", ctypes.c_float).reshape(shapes[0])
]
wavelets = [np.fromfile("data_2_wavelet.bin", ctypes.c_float)]
wavelet_idxs = [np.fromfile("data_2_wavelet_idxs.bin", ctypes.c_int)]
initial_factor = 1.0
n_its = 1250
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
    wavelet_idxs=wavelet_idxs,
    wavelets=wavelets,
)

data[0].tofile("out/data_2_deblended_data.bin")
