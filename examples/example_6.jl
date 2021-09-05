using MPI

include("../wrappers/julia/agdeblend.jl")
using .AGDeblend

MPI.Init()

comm = MPI.COMM_WORLD
comm_rank = MPI.Comm_rank(comm)

if comm_rank == 0
  n_traces = [80]
else
  n_traces = [64]
end
n_times = 512
volumes = [0]
n_dims = [2]
window_shapes = [[256, 32]]
coords = [[comm_rank]]
shottimes = [zeros(Clong, n_traces[1])]
open("data_1_shottimes.bin", "r") do f
  seek(f, comm_rank * 64 * sizeof(Clong))
  read!(f, shottimes[1])
end
channels = [zeros(Cint, n_traces[1])]
open("data_1_channels.bin", "r") do f
  seek(f, comm_rank * 64 * sizeof(Cint))
  read!(f, channels[1])
end
trace_types = [zeros(Cint, n_traces[1])]
open("data_1_trace_types.bin", "r") do f
  seek(f, comm_rank * 64 * sizeof(Cint))
  read!(f, trace_types[1])
end
data = [zeros(Cfloat, (n_times, n_traces[1]))]
open("data_1_blended_data.bin", "r") do f
  seek(f, comm_rank * 64 * n_times * sizeof(Cfloat))
  read!(f, data[1])
end
initial_factor = 1.0
n_its = 1000
print_freq = -1

deblend(volumes, window_shapes, coords, shottimes, channels, trace_types,
        initial_factor, n_its, print_freq, data, comm=comm)

n_overlap = div(window_shapes[1][2], 2)

if comm_rank == 0
  buf = MPI.recv(1, 0, comm)[1]
  data[1][:, n_traces[1] - n_overlap + 1:end,] += buf
  write("out/data_1_deblended_data.bin", data[1])
  MPI.Barrier(comm)
else
  MPI.send(data[1][:, 1:n_overlap], 0, 0, comm)
  MPI.Barrier(comm)
  open("out/data_1_deblended_data.bin", "a") do f
    write(f, data[1][:, n_overlap + 1:end])
  end
end

MPI.Finalize()
