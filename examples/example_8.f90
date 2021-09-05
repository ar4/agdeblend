program example_8

use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float
use :: agdeblend_m, only: deblend, deblend_patch, volume

implicit none

type(deblend_patch), dimension(:), allocatable :: patches
type(volume), dimension(:), allocatable :: volumes
integer, parameter :: n_traces = 128
integer, parameter :: n_times(2) = (/ 512, 768 /)
integer, parameter :: n_window_traces = 64
integer, parameter :: n_window_times = 256
real, parameter :: initial_factor = 1.0
integer, parameter :: n_its = 2500
integer, parameter :: print_freq = -1
integer :: fid
integer :: volume_idx
character(len=256) :: filename
integer :: ierr

! Set patches
allocate(patches(2))
allocate(volumes(2))
do volume_idx = 1, 2
  patches(volume_idx)%volume_idx = volume_idx - 1
  allocate(patches(volume_idx)%coords(1))
  patches(volume_idx)%coords(1) = 0
  allocate(patches(volume_idx)%data_shape(2))
  patches(volume_idx)%data_shape(1) = n_times(volume_idx)
  patches(volume_idx)%data_shape(2) = n_traces
  allocate(patches(volume_idx)%shottimes(n_traces))
  allocate(patches(volume_idx)%channels(n_traces))
  allocate(patches(volume_idx)%trace_types(n_traces))
  allocate(patches(volume_idx)%values(n_traces * n_times(volume_idx)))
  volumes(volume_idx)%n_dims = 2
  allocate(volumes(volume_idx)%window_shape(2))
  volumes(volume_idx)%window_shape(1) = n_window_times
  volumes(volume_idx)%window_shape(2) = n_window_traces

  write (filename, '(a, i0, a)') "data_3_shottimes_", volume_idx - 1,   &
    ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="read", status="old")
  read(fid) patches(volume_idx)%shottimes
  close(fid)

  write (filename, '(a, i0, a)') "data_3_channels_", volume_idx - 1,    &
    ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="read", status="old")
  read(fid) patches(volume_idx)%channels
  close(fid)

  write (filename, '(a, i0, a)') "data_3_trace_types_", volume_idx - 1, &
    ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="read", status="old")
  read(fid) patches(volume_idx)%trace_types
  close(fid)

  write (filename, '(a, i0, a)') "data_3_blended_data_", volume_idx - 1,&
    ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="read", status="old")
  read(fid) patches(volume_idx)%values
  close(fid)
end do

call deblend(patches, volumes, initial_factor, n_its, print_freq, ierr)

if (ierr /= 0) then
  stop 1
end if

do volume_idx = 1, 2
  write (filename, '(a, i0, a)') "out/data_3_deblended_data_",          &
    volume_idx - 1, ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="write", status="replace")
  write(fid) patches(volume_idx)%values
  close(fid)

  deallocate(patches(volume_idx)%coords)
  deallocate(patches(volume_idx)%data_shape)
  deallocate(patches(volume_idx)%shottimes)
  deallocate(patches(volume_idx)%channels)
  deallocate(patches(volume_idx)%trace_types)
  deallocate(patches(volume_idx)%values)
  deallocate(volumes(volume_idx)%window_shape)
end do
deallocate(patches)
deallocate(volumes)

end program example_8
