import ctypes
import numpy as np
from agdeblend import blend, AGDBlendMean

n_traces = [128]
n_times = [512]
n_times_out = [768]
shottimes = [np.fromfile("data_1_shottimes.bin", ctypes.c_long)]
shottimes_out = [np.fromfile("data_1_shottimes_out.bin", ctypes.c_long)]
channels = [np.fromfile("data_1_channels.bin", ctypes.c_int)]
trace_types = [np.fromfile("data_1_trace_types.bin", ctypes.c_int)]
data = [
    np.fromfile("data_1_blended_data.bin", ctypes.c_float).reshape(
        n_traces[0], n_times[0]
    )
]
blend_mode = AGDBlendMean
taper_length = 16

blended_data = blend(
    shottimes,
    channels,
    trace_types,
    data,
    blend_mode,
    taper_length=taper_length,
    n_times_out=n_times_out,
    shottimes_out=shottimes_out,
)

blended_data[0].tofile("out/data_1_blended_data_2.bin")
