using MPI

include("../wrappers/julia/agdeblend.jl")
using .AGDeblend

MPI.Init()

comm = MPI.COMM_WORLD
comm_rank = MPI.Comm_rank(comm)

n_traces = [64]
n_times = [512]
shottimes = [zeros(Clong, n_traces[1])]
open("data_1_shottimes.bin", "r") do f
  seek(f, comm_rank * n_traces[1] * sizeof(Clong))
  read!(f, shottimes[1])
end
channels = [zeros(Cint, n_traces[1])]
open("data_1_channels.bin", "r") do f
  seek(f, comm_rank * n_traces[1] * sizeof(Cint))
  read!(f, channels[1])
end
trace_types = [zeros(Cint, n_traces[1])]
open("data_1_trace_types.bin", "r") do f
  seek(f, comm_rank * n_traces[1] * sizeof(Cint))
  read!(f, trace_types[1])
end
data = [zeros(Cfloat, (n_times[1], n_traces[1]))]
open("data_1_true_data.bin", "r") do f
  seek(f, comm_rank * n_traces[1] * n_times[1] * sizeof(Cfloat))
  read!(f, data[1])
end
blend_mode = AGDBlendSum

blended_data = blend(shottimes, channels, trace_types, data, blend_mode,
                     comm=comm)
if comm_rank == 0
  write("out/data_1_blended_data.bin", blended_data[1])
  MPI.Barrier(comm)
else
  MPI.Barrier(comm)
  open("out/data_1_blended_data.bin", "a") do f
    write(f, blended_data[1])
  end
end

MPI.Finalize()
