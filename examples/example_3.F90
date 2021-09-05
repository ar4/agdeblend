program example_3

use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float, c_sizeof
#ifdef AGD_F08
use :: mpi_f08, only: MPI_Init, MPI_Comm_rank, MPI_COMM_WORLD,          &
  MPI_Barrier, MPI_Finalize
#else
use :: mpi, only: MPI_Init, MPI_Comm_rank, MPI_COMM_WORLD,              &
  MPI_Barrier, MPI_Finalize
#endif /* AGD_F08 */
use :: agdeblend_m, only: blend, blend_patch, AGDBlendSum

implicit none

type(blend_patch), dimension(:), allocatable :: patches
integer, parameter :: n_traces = 64
integer, parameter :: n_times = 512
integer, parameter :: blend_mode = AGDBlendSum
integer :: fid
integer :: comm_rank
integer :: ierr

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, comm_rank, ierr)

! Set patches
allocate(patches(1))
patches(1)%n_traces = n_traces
patches(1)%n_times = n_times
allocate(patches(1)%shottimes(n_traces))
allocate(patches(1)%channels(n_traces))
allocate(patches(1)%trace_types(n_traces))
allocate(patches(1)%values(n_traces * n_times))

open(newunit=fid, file="data_1_shottimes.bin", form="unformatted",      &
     access="stream", action="read", status="old")
read(fid,                                                               &
     pos=1 + comm_rank * c_sizeof(patches(1)%shottimes(1)) * n_traces)  &
  patches(1)%shottimes
close(fid)
open(newunit=fid, file="data_1_channels.bin", form="unformatted",       &
     access="stream", action="read", status="old")
read(fid, pos=1 + c_sizeof(patches(1)%channels(1)) * n_traces)          &
  patches(1)%channels
close(fid)
open(newunit=fid, file="data_1_trace_types.bin", form="unformatted",    &
     access="stream", action="read", status="old")
read(fid, pos=1 + c_sizeof(patches(1)%trace_types(1)) * n_traces)       &
  patches(1)%trace_types
close(fid)
open(newunit=fid, file="data_1_true_data.bin", form="unformatted",      &
     access="stream", action="read", status="old")
read(fid, pos=1 + c_sizeof(patches(1)%values(1)) * n_traces * n_times)  &
  patches(1)%values
close(fid)

call blend(patches, blend_mode, MPI_COMM_WORLD, ierr)

if (ierr /= 0) then
  stop 1
end if

if (comm_rank == 0) then
  open(newunit=fid, file="out/data_1_blended_data.bin",                 &
       form="unformatted", access="stream", action="write",             &
       status="replace")
  write(fid) patches(1)%values
  close(fid)
  call MPI_Barrier(MPI_COMM_WORLD, ierr);
else
  call MPI_Barrier(MPI_COMM_WORLD, ierr);
  open(newunit=fid, file="out/data_1_blended_data.bin",                 &
       form="unformatted", access="stream", action="readwrite",         &
       status="old", position="append")
  write(fid) patches(1)%values
  close(fid)
end if

deallocate(patches(1)%shottimes)
deallocate(patches(1)%channels)
deallocate(patches(1)%trace_types)
deallocate(patches(1)%values)
deallocate(patches)

call MPI_Finalize(ierr)

end program example_3
