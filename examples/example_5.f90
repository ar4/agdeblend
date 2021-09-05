program example_5

use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float, c_sizeof
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
integer :: ierr

! Set patches
allocate(patches(2))
patches(1)%volume_idx = 0
patches(2)%volume_idx = 0
allocate(patches(1)%coords(1))
allocate(patches(2)%coords(1))
patches(1)%coords(1) = 0
patches(2)%coords(1) = 1
allocate(patches(1)%data_shape(2))
allocate(patches(2)%data_shape(2))
patches(1)%data_shape(1) = n_times
patches(2)%data_shape(1) = n_times
patches(1)%data_shape(2) = n_traces0
patches(2)%data_shape(2) = n_traces1
allocate(patches(1)%shottimes(n_traces0))
allocate(patches(2)%shottimes(n_traces1))
allocate(patches(1)%channels(n_traces0))
allocate(patches(2)%channels(n_traces1))
allocate(patches(1)%trace_types(n_traces0))
allocate(patches(2)%trace_types(n_traces1))
allocate(patches(1)%values(n_traces0 * n_times))
allocate(patches(2)%values(n_traces1 * n_times))
allocate(volumes(1))
volumes(1)%n_dims = 2
allocate(volumes(1)%window_shape(2))
volumes(1)%window_shape(1) = n_window_times
volumes(1)%window_shape(2) = n_window_traces

open(newunit=fid, file="data_1_shottimes.bin", form="unformatted",      &
     access="stream", action="read", status="old")
read(fid) patches(1)%shottimes
read(fid, pos=1 + c_sizeof(patches(2)%shottimes(1)) * 64)               &
  patches(2)%shottimes
close(fid)

open(newunit=fid, file="data_1_channels.bin", form="unformatted",       &
     access="stream", action="read", status="old")
read(fid) patches(1)%channels
read(fid, pos=1 + c_sizeof(patches(2)%channels(1)) * 64)                &
  patches(2)%channels
close(fid)

open(newunit=fid, file="data_1_trace_types.bin", form="unformatted",    &
     access="stream", action="read", status="old")
read(fid) patches(1)%trace_types
read(fid, pos=1 + c_sizeof(patches(2)%trace_types(1)) * 64)             &
  patches(2)%trace_types
close(fid)

open(newunit=fid, file="data_1_blended_data.bin", form="unformatted",   &
     access="stream", action="read", status="old")
read(fid) patches(1)%values
read(fid, pos=1 + c_sizeof(patches(2)%values(1)) * 64 * n_times)        &
  patches(2)%values
close(fid)

call deblend(patches, volumes, initial_factor, n_its, print_freq, ierr)

if (ierr /= 0) then
  stop 1
end if

! Add overlap between patches 1 and 2 to 1
do trace_idx = 1, n_overlap
  trace0_idx = n_traces0 - n_overlap + trace_idx
  do time_idx = 1, n_times
    patches(1)%values(trace0_idx * n_times + time_idx) =                &
      patches(1)%values(trace0_idx * n_times + time_idx) +              &
      patches(2)%values(trace_idx * n_times + time_idx)
  end do
end do

open(newunit=fid, file="out/data_1_deblended_data.bin",                 &
     form="unformatted", access="stream", action="write",               &
     status="replace")
write(fid) patches(1)%values
write(fid) patches(2)%values(n_overlap * n_times + 1:)
close(fid)

deallocate(patches(1)%coords)
deallocate(patches(2)%coords)
deallocate(patches(1)%data_shape)
deallocate(patches(2)%data_shape)
deallocate(patches(1)%shottimes)
deallocate(patches(2)%shottimes)
deallocate(patches(1)%channels)
deallocate(patches(2)%channels)
deallocate(patches(1)%trace_types)
deallocate(patches(2)%trace_types)
deallocate(patches(1)%values)
deallocate(patches(2)%values)
deallocate(patches)
deallocate(volumes(1)%window_shape)
deallocate(volumes)

end program example_5
