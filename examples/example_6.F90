program example_6

use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float, c_sizeof
#ifdef AGD_F08
use :: mpi_f08, only: MPI_Init, MPI_Comm_rank, MPI_COMM_WORLD,          &
  MPI_Barrier, MPI_Finalize, MPI_Send, MPI_Recv, MPI_FLOAT,             &
  MPI_STATUS_IGNORE
#else
use :: mpi, only: MPI_Init, MPI_Comm_rank, MPI_COMM_WORLD,              &
  MPI_Barrier, MPI_Finalize, MPI_Send, MPI_Recv, MPI_FLOAT,             &
  MPI_STATUS_IGNORE
#endif /* AGD_F08 */

use :: agdeblend_m, only: deblend, deblend_patch, volume

implicit none

type(deblend_patch), dimension(:), allocatable :: patches
type(volume), dimension(:), allocatable :: volumes
integer, parameter :: n_traces0 = 80
integer, parameter :: n_traces1 = 64
integer, parameter :: n_times = 512
integer, parameter :: n_window_traces = 32
integer, parameter :: n_window_times = 256
real, parameter :: initial_factor = 1.0
integer, parameter :: n_its = 1000
integer, parameter :: print_freq = -1
integer, parameter :: n_overlap = n_window_traces / 2
integer :: fid
integer :: trace_idx
integer :: trace0_idx
integer :: time_idx
integer :: n_traces
real, allocatable :: buffer(:)
integer :: comm_rank
integer :: ierr

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, comm_rank, ierr)

if (comm_rank == 0) then
  n_traces = n_traces0
else
  n_traces = n_traces1
end if

! Set patches
allocate(patches(1))
patches(1)%volume_idx = 0
allocate(patches(1)%coords(1))
patches(1)%coords(1) = comm_rank
allocate(patches(1)%data_shape(2))
patches(1)%data_shape(1) = n_times
patches(1)%data_shape(2) = n_traces
allocate(patches(1)%shottimes(n_traces))
allocate(patches(1)%channels(n_traces))
allocate(patches(1)%trace_types(n_traces))
allocate(patches(1)%values(n_traces * n_times))
allocate(volumes(1))
volumes(1)%n_dims = 2
allocate(volumes(1)%window_shape(2))
volumes(1)%window_shape(1) = n_window_times
volumes(1)%window_shape(2) = n_window_traces

open(newunit=fid, file="data_1_shottimes.bin", form="unformatted",      &
     access="stream", action="read", status="old")
read(fid, pos=1 + comm_rank * c_sizeof(patches(1)%shottimes(1)) * 64)   &
  patches(1)%shottimes
close(fid)

open(newunit=fid, file="data_1_channels.bin", form="unformatted",       &
     access="stream", action="read", status="old")
read(fid, pos=1 + comm_rank * c_sizeof(patches(1)%channels(1)) * 64)    &
  patches(1)%channels
close(fid)

open(newunit=fid, file="data_1_trace_types.bin", form="unformatted",    &
     access="stream", action="read", status="old")
read(fid, pos=1 + comm_rank * c_sizeof(patches(1)%trace_types(1)) * 64) &
  patches(1)%trace_types
close(fid)

open(newunit=fid, file="data_1_blended_data.bin", form="unformatted",   &
     access="stream", action="read", status="old")
read(fid,                                                               &
     pos=1 + comm_rank * c_sizeof(patches(2)%values(1)) * 64 * n_times) &
  patches(1)%values
close(fid)

call deblend(patches, volumes, initial_factor, n_its, print_freq,       &
             MPI_COMM_WORLD, ierr)

if (ierr /= 0) then
  stop 1
end if

! Add overlap between patches 1 and 2 to 1 and write to file
if (comm_rank == 0) then
  allocate(buffer(n_overlap * n_times))
  call MPI_Recv(buffer, n_overlap * n_times, MPI_FLOAT, 1, 0,           &
                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  do trace_idx = 1, n_overlap
    trace0_idx = n_traces0 - n_overlap + trace_idx
    do time_idx = 1, n_times
      patches(1)%values(trace0_idx * n_times + time_idx) =              &
        patches(1)%values(trace0_idx * n_times + time_idx) +            &
        buffer(trace_idx * n_times + time_idx)
    end do
  end do
  open(newunit=fid, file="out/data_1_deblended_data.bin",               &
       form="unformatted", access="stream", action="write",             &
       status="replace")
  write(fid) patches(1)%values
  close(fid)
  deallocate(buffer)
  call MPI_Barrier(MPI_COMM_WORLD, ierr);
else
  call MPI_Send(patches(1)%values, n_overlap * n_times, MPI_FLOAT, 0, 0,&
                MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr);
  open(newunit=fid, file="out/data_1_deblended_data.bin",               &
       form="unformatted", access="stream", action="readwrite",         &
       status="old", position="append")
  write(fid) patches(1)%values(n_overlap * n_times + 1:)
  close(fid)
end if

deallocate(patches(1)%coords)
deallocate(patches(1)%data_shape)
deallocate(patches(1)%shottimes)
deallocate(patches(1)%channels)
deallocate(patches(1)%trace_types)
deallocate(patches(1)%values)
deallocate(patches)
deallocate(volumes(1)%window_shape)
deallocate(volumes)

call MPI_Finalize(ierr)

end program example_6
