program example_1

use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float
use :: agdeblend_m, only: blend, blend_patch, AGDBlendSum

implicit none

type(blend_patch), dimension(:), allocatable :: patches
integer, parameter :: n_traces = 128
integer, parameter :: n_times = 512
integer, parameter :: blend_mode = AGDBlendSum
integer :: fid
integer :: ierr

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
read(fid) patches(1)%shottimes
close(fid)
open(newunit=fid, file="data_1_channels.bin", form="unformatted",       &
     access="stream", action="read", status="old")
read(fid) patches(1)%channels
close(fid)
open(newunit=fid, file="data_1_trace_types.bin", form="unformatted",    &
     access="stream", action="read", status="old")
read(fid) patches(1)%trace_types
close(fid)
open(newunit=fid, file="data_1_true_data.bin", form="unformatted",      &
     access="stream", action="read", status="old")
read(fid) patches(1)%values
close(fid)

call blend(patches, blend_mode, ierr)

if (ierr /= 0) then
  stop 1
end if

open(newunit=fid, file="out/data_1_blended_data.bin",                   &
     form="unformatted", access="stream", action="write",               &
     status="replace")
write(fid) patches(1)%values
close(fid)

deallocate(patches(1)%shottimes)
deallocate(patches(1)%channels)
deallocate(patches(1)%trace_types)
deallocate(patches(1)%values)
deallocate(patches)

end program example_1
