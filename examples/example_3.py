import ctypes
import numpy as np
from mpi4py import MPI
from agdeblend import blend, AGDBlendSum

comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()

n_traces = [64]
n_times = [512]
shottimes = [
    np.fromfile(
        "data_1_shottimes.bin",
        ctypes.c_long,
        offset=comm_rank * n_traces[0] * ctypes.sizeof(ctypes.c_long),
        count=n_traces[0],
    )
]
channels = [
    np.fromfile(
        "data_1_channels.bin",
        ctypes.c_int,
        offset=comm_rank * n_traces[0] * ctypes.sizeof(ctypes.c_int),
        count=n_traces[0],
    )
]
trace_types = [
    np.fromfile(
        "data_1_trace_types.bin",
        ctypes.c_int,
        offset=comm_rank * n_traces[0] * ctypes.sizeof(ctypes.c_int),
        count=n_traces[0],
    )
]
data = [
    np.fromfile(
        "data_1_true_data.bin",
        ctypes.c_float,
        offset=comm_rank
        * n_traces[0]
        * n_times[0]
        * ctypes.sizeof(ctypes.c_float),
        count=n_traces[0] * n_times[0],
    ).reshape(n_traces[0], n_times[0])
]
blend_mode = AGDBlendSum

blended_data = blend(
    shottimes, channels, trace_types, data, blend_mode, comm=comm
)

if comm_rank == 0:
    blended_data[0].tofile("out/data_1_blended_data.bin")
    comm.Barrier()
else:
    comm.Barrier()
    with open("out/data_1_blended_data.bin", mode="a") as fid:
        blended_data[0].tofile(fid)

MPI.Finalize()
