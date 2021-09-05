program example_7

use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float
use :: agdeblend_m, only: deblend, deblend_patch, volume, wavelet

implicit none

type(deblend_patch), dimension(:), allocatable :: patches
type(volume), dimension(:), allocatable :: volumes
type(wavelet), dimension(:), allocatable :: wavelets
integer, parameter :: n_traces = 64
integer, parameter :: n_times = 1024
integer, parameter :: n_window_traces = 32
integer, parameter :: n_window_times = 256
real, parameter :: initial_factor = 1.0
integer, parameter :: n_its = 1250
integer, parameter :: print_freq = -1
integer :: fid
integer :: ierr

! Set patches
allocate(patches(1))
patches(1)%volume_idx = 0
allocate(patches(1)%coords(1))
patches(1)%coords(1) = 0
allocate(patches(1)%data_shape(2))
patches(1)%data_shape(1) = n_times
patches(1)%data_shape(2) = n_traces
allocate(patches(1)%shottimes(n_traces))
allocate(patches(1)%channels(n_traces))
allocate(patches(1)%trace_types(n_traces))
allocate(patches(1)%wavelet_idxs(n_traces))
allocate(patches(1)%values(n_traces * n_times))

! Set volumes
allocate(volumes(1))
volumes(1)%n_dims = 2
allocate(volumes(1)%window_shape(2))
volumes(1)%window_shape(1) = n_window_times
volumes(1)%window_shape(2) = n_window_traces

! Set wavelets
allocate(wavelets(1))
allocate(wavelets(1)%values(513))

open(newunit=fid, file="data_2_shottimes.bin", form="unformatted",      &
     access="stream", action="read", status="old")
read(fid) patches(1)%shottimes
close(fid)
open(newunit=fid, file="data_2_channels.bin", form="unformatted",       &
     access="stream", action="read", status="old")
read(fid) patches(1)%channels
close(fid)
open(newunit=fid, file="data_2_trace_types.bin", form="unformatted",    &
     access="stream", action="read", status="old")
read(fid) patches(1)%trace_types
close(fid)
open(newunit=fid, file="data_2_blended_data.bin", form="unformatted",   &
     access="stream", action="read", status="old")
read(fid) patches(1)%values
close(fid)
open(newunit=fid, file="data_2_wavelet.bin", form="unformatted",        &
     access="stream", action="read", status="old")
read(fid) wavelets(1)%values
close(fid)
open(newunit=fid, file="data_2_wavelet_idxs.bin", form="unformatted",   &
     access="stream", action="read", status="old")
read(fid) patches(1)%wavelet_idxs
close(fid)

call deblend(patches, volumes, initial_factor, n_its, print_freq, ierr, &
             wavelets)

if (ierr /= 0) then
  stop 1
end if

open(newunit=fid, file="out/data_2_deblended_data.bin",                 &
     form="unformatted", access="stream", action="write",               &
     status="replace")
write(fid) patches(1)%values
close(fid)

deallocate(patches(1)%coords)
deallocate(patches(1)%data_shape)
deallocate(patches(1)%shottimes)
deallocate(patches(1)%channels)
deallocate(patches(1)%trace_types)
deallocate(patches(1)%wavelet_idxs)
deallocate(patches(1)%values)
deallocate(patches)
deallocate(volumes(1)%window_shape)
deallocate(volumes)
deallocate(wavelets(1)%values)
deallocate(wavelets)

end program example_7
