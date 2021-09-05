import ctypes
from mpi4py import MPI
import numpy as np
from agdeblend import deblend

comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()

if comm_rank == 0:
    n_traces = [80]
else:
    n_traces = [64]
n_times = 512
volumes = [0]
n_dims = [2]
window_shapes = [[32, 256]]
coords = [[comm_rank]]
shapes = [[n_traces[0], n_times]]
shottimes = [
    np.fromfile(
        "data_1_shottimes.bin",
        ctypes.c_long,
        offset=comm_rank * 64 * ctypes.sizeof(ctypes.c_long),
        count=n_traces[0],
    )
]
channels = [
    np.fromfile(
        "data_1_channels.bin",
        ctypes.c_int,
        offset=comm_rank * 64 * ctypes.sizeof(ctypes.c_int),
        count=n_traces[0],
    )
]
trace_types = [
    np.fromfile(
        "data_1_trace_types.bin",
        ctypes.c_int,
        offset=comm_rank * 64 * ctypes.sizeof(ctypes.c_int),
        count=n_traces[0],
    )
]
data = [
    np.fromfile(
        "data_1_blended_data.bin",
        ctypes.c_float,
        offset=comm_rank * 64 * n_times * ctypes.sizeof(ctypes.c_float),
        count=n_traces[0] * n_times,
    ).reshape(shapes[0])
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
    comm=comm,
)

n_overlap = window_shapes[0][0] // 2

if comm_rank == 0:
    buf = comm.recv(source=1)
    data[0][n_traces[0] - n_overlap :] += buf
    with open("out/data_1_deblended_data.bin", mode="w") as fid:
        data[0].tofile(fid)
    comm.Barrier()
else:
    comm.send(data[0][:n_overlap], dest=0)
    comm.Barrier()
    with open("out/data_1_deblended_data.bin", mode="a") as fid:
        data[0][n_overlap:].tofile(fid)
