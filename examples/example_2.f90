program example_2

use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float
use :: agdeblend_m, only: blend, blend_patch, AGDBlendMean

implicit none

type(blend_patch), dimension(:), allocatable :: patches
type(blend_patch), dimension(:), allocatable :: patches_out
integer, parameter :: n_traces = 128
integer, parameter :: n_times = 512
integer, parameter :: n_times_out = 768
integer, parameter :: taper_length = 16
integer, parameter :: blend_mode = AGDBlendMean
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

allocate(patches_out(1))
patches_out(1)%n_traces = n_traces
patches_out(1)%n_times = n_times_out
allocate(patches_out(1)%shottimes(n_traces))
allocate(patches_out(1)%channels(n_traces))
allocate(patches_out(1)%trace_types(n_traces))
allocate(patches_out(1)%values(n_traces * n_times_out))

open(newunit=fid, file="data_1_shottimes.bin", form="unformatted",      &
     access="stream", action="read", status="old")
read(fid) patches(1)%shottimes
close(fid)
open(newunit=fid, file="data_1_shottimes_out.bin", form="unformatted",  &
     access="stream", action="read", status="old")
read(fid) patches_out(1)%shottimes
close(fid)
open(newunit=fid, file="data_1_channels.bin", form="unformatted",       &
     access="stream", action="read", status="old")
read(fid) patches(1)%channels
close(fid)
open(newunit=fid, file="data_1_channels.bin", form="unformatted",       &
     access="stream", action="read", status="old")
read(fid) patches_out(1)%channels
close(fid)
open(newunit=fid, file="data_1_trace_types.bin", form="unformatted",    &
     access="stream", action="read", status="old")
read(fid) patches(1)%trace_types
close(fid)
open(newunit=fid, file="data_1_trace_types.bin", form="unformatted",    &
     access="stream", action="read", status="old")
read(fid) patches_out(1)%trace_types
close(fid)
open(newunit=fid, file="data_1_blended_data.bin", form="unformatted",   &
     access="stream", action="read", status="old")
read(fid) patches(1)%values
close(fid)

call blend(patches, blend_mode, ierr, taper_length=taper_length,        &
           patches_out=patches_out)

if (ierr /= 0) then
  stop 1
end if

open(newunit=fid, file="out/data_1_blended_data_2.bin",                 &
     form="unformatted", access="stream", action="write",               &
     status="replace")
write(fid) patches_out(1)%values
close(fid)

deallocate(patches(1)%shottimes)
deallocate(patches(1)%channels)
deallocate(patches(1)%trace_types)
deallocate(patches(1)%values)
deallocate(patches)

deallocate(patches_out(1)%shottimes)
deallocate(patches_out(1)%channels)
deallocate(patches_out(1)%trace_types)
deallocate(patches_out(1)%values)
deallocate(patches_out)

end program example_2
