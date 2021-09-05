program example_9

use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float
use :: agdeblend_m, only: deblend, deblend_patch, volume, AGDMissing

implicit none

type(deblend_patch), dimension(:), allocatable :: patches
type(volume), dimension(:), allocatable :: volumes
integer, parameter :: n_times = 512
integer, parameter :: n_patches = 7
integer, parameter :: n_window_traces_x = 16
integer, parameter :: n_window_traces_y = 16
integer, parameter :: n_window_times = 256
real, parameter :: initial_factor = 1.0
integer, parameter :: n_its = 1000
integer, parameter :: print_freq = -1
integer :: fid
integer :: patch_idx
integer :: file_idx
integer :: ix
integer :: iy
character(len=256) :: filename
integer :: ierr
integer, parameter :: space_shape_all(2) = (/ 40, 32 /)
integer, parameter :: file_ranges(2, 2, 2) =                            &
  reshape(source = (/ 1, 24, 1, 16, 17, 40, 17, 32 /),                  &
          shape = shape(file_ranges))
integer(kind=c_long), allocatable :: shottimes_all(:, :)
integer(kind=c_int), allocatable :: channels_all(:, :)
integer(kind=c_int), allocatable :: trace_types_all(:, :)
real(kind=c_float), allocatable :: data_all(:, :, :)
integer, parameter :: patch_ranges_y(2, 3) =                            &
  reshape(source = (/ 1, 16, 9, 32, 25, 40 /),                          &
          shape = shape(patch_ranges_y))
integer, parameter :: patch_ranges_x(2, 3) =                            &
  reshape(source = (/ 1, 16, 9, 24, 17, 32 /),                          &
          shape = shape(patch_ranges_x))
integer :: patch_ranges(2, 2, 7)
integer :: n_traces_y
integer :: n_traces_x

allocate(shottimes_all(space_shape_all(1), space_shape_all(2)))
allocate(channels_all(space_shape_all(1), space_shape_all(2)))
allocate(trace_types_all(space_shape_all(1), space_shape_all(2)))
allocate(data_all(n_times, space_shape_all(1), space_shape_all(2)))

trace_types_all = AGDMissing

! Read input files
do file_idx = 1, 2
  write (filename, '(a, i0, a)') "data_4_shottimes_", file_idx - 1,     &
    ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="read", status="old")
  read(fid) shottimes_all(file_ranges(1, 1, file_idx):                  &
                          file_ranges(2, 1, file_idx),                  &
                          file_ranges(1, 2, file_idx):                  &
                          file_ranges(2, 2, file_idx))
  close(fid)

  write (filename, '(a, i0, a)') "data_4_channels_", file_idx - 1,      &
    ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="read", status="old")
  read(fid) channels_all(file_ranges(1, 1, file_idx):                   &
                         file_ranges(2, 1, file_idx),                   &
                         file_ranges(1, 2, file_idx):                   &
                         file_ranges(2, 2, file_idx))
  close(fid)

  write (filename, '(a, i0, a)') "data_4_trace_types_", file_idx - 1,   &
    ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="read", status="old")
  read(fid) trace_types_all(file_ranges(1, 1, file_idx):                &
                            file_ranges(2, 1, file_idx),                &
                            file_ranges(1, 2, file_idx):                &
                            file_ranges(2, 2, file_idx))
  close(fid)

  write (filename, '(a, i0, a)') "data_4_blended_data_", file_idx - 1,  &
    ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="read", status="old")
  read(fid) data_all(:,                                                 &
                     file_ranges(1, 1, file_idx):                       &
                     file_ranges(2, 1, file_idx),                       &
                     file_ranges(1, 2, file_idx):                       &
                     file_ranges(2, 2, file_idx))
  close(fid)
end do

! Set patches
allocate(patches(n_patches))
patch_idx = 1
do ix = 1, 3
  do iy = 1, 3
    if ((ix == 3 .and. iy == 1) .or. (ix == 1 .and. iy == 3)) then
      cycle
    end if
    patch_ranges(:, 1, patch_idx) = patch_ranges_y(:, iy)
    patch_ranges(:, 2, patch_idx) = patch_ranges_x(:, ix)
    n_traces_y = patch_ranges_y(2, iy) - patch_ranges_y(1, iy) + 1
    n_traces_x = patch_ranges_x(2, ix) - patch_ranges_x(1, ix) + 1
    patches(patch_idx)%volume_idx = 0
    allocate(patches(patch_idx)%coords(2))
    patches(patch_idx)%coords(1) = iy
    patches(patch_idx)%coords(2) = ix
    allocate(patches(patch_idx)%data_shape(3))
    patches(patch_idx)%data_shape(1) = n_times
    patches(patch_idx)%data_shape(2) = n_traces_y
    patches(patch_idx)%data_shape(3) = n_traces_x
    allocate(patches(patch_idx)%shottimes(n_traces_y * n_traces_x))
    allocate(patches(patch_idx)%channels(n_traces_y * n_traces_x))
    allocate(patches(patch_idx)%trace_types(n_traces_y * n_traces_x))
    allocate(patches(patch_idx)%values                                  &
      (n_times * n_traces_y * n_traces_x))
    patches(patch_idx)%shottimes(:) =                                   &
      reshape(source = shottimes_all(patch_ranges(1, 1, patch_idx) :    &
                                     patch_ranges(2, 1, patch_idx),     &
                                     patch_ranges(1, 2, patch_idx) :    &
                                     patch_ranges(2, 2, patch_idx)),    &
              shape = (/ n_traces_y * n_traces_x /))
    patches(patch_idx)%channels(:) =                                    &
      reshape(source = channels_all(patch_ranges(1, 1, patch_idx) :     &
                                    patch_ranges(2, 1, patch_idx),      &
                                    patch_ranges(1, 2, patch_idx) :     &
                                    patch_ranges(2, 2, patch_idx)),     &
              shape = (/ n_traces_y * n_traces_x /))
    patches(patch_idx)%trace_types(:) =                                 &
      reshape(source = trace_types_all(patch_ranges(1, 1, patch_idx) :  &
                                       patch_ranges(2, 1, patch_idx),   &
                                       patch_ranges(1, 2, patch_idx) :  &
                                       patch_ranges(2, 2, patch_idx)),  &
              shape = (/ n_traces_y * n_traces_x /))
    patches(patch_idx)%values(:) =                                      &
      reshape(source = data_all(:,                                      &
                                patch_ranges(1, 1, patch_idx) :         &
                                patch_ranges(2, 1, patch_idx),          &
                                patch_ranges(1, 2, patch_idx) :         &
                                patch_ranges(2, 2, patch_idx)),         &
              shape = (/ n_times * n_traces_y * n_traces_x /))
    patch_idx = patch_idx + 1
  end do
end do

deallocate(shottimes_all)
deallocate(channels_all)
deallocate(trace_types_all)
deallocate(data_all)

allocate(volumes(1))
volumes(1)%n_dims = 3
allocate(volumes(1)%window_shape(3))
volumes(1)%window_shape(1) = n_window_times
volumes(1)%window_shape(2) = n_window_traces_y
volumes(1)%window_shape(3) = n_window_traces_x

! Deblend
call deblend(patches, volumes, initial_factor, n_its, print_freq, ierr)

if (ierr /= 0) then
  stop 1
end if

! Sum patches
allocate(data_all(n_times, space_shape_all(1), space_shape_all(2)))
data_all = 0

do patch_idx = 1, n_patches
  n_traces_y = patch_ranges(2, 1, patch_idx) -                          &
               patch_ranges(1, 1, patch_idx) + 1
  n_traces_x = patch_ranges(2, 2, patch_idx) -                          &
               patch_ranges(1, 2, patch_idx) + 1
  data_all(:,                                                           &
           patch_ranges(1, 1, patch_idx):patch_ranges(2, 1, patch_idx), &
           patch_ranges(1, 2, patch_idx):patch_ranges(2, 2, patch_idx)) &
    = data_all(:,                                                       &
           patch_ranges(1, 1, patch_idx):patch_ranges(2, 1, patch_idx), &
           patch_ranges(1, 2, patch_idx):patch_ranges(2, 2, patch_idx)) &
    + reshape(source = patches(patch_idx)%values,                       &
              shape = (/ n_times, n_traces_y, n_traces_x /))
end do

! Extract file ranges and save to disk
do file_idx = 1, 2
  write (filename, '(a, i0, a)') "out/data_4_deblended_data_",          &
         file_idx - 1, ".bin"
  open(newunit=fid, file=trim(filename), form="unformatted",            &
       access="stream", action="write", status="replace")
  write(fid) data_all(:,                                                &
                      file_ranges(1, 1, file_idx):                      &
                      file_ranges(2, 1, file_idx),                      &
                      file_ranges(1, 2, file_idx):                      &
                      file_ranges(2, 2, file_idx))
  close(fid)
end do
deallocate(data_all)

do patch_idx = 1, n_patches
  deallocate(patches(patch_idx)%coords)
  deallocate(patches(patch_idx)%data_shape)
  deallocate(patches(patch_idx)%shottimes)
  deallocate(patches(patch_idx)%channels)
  deallocate(patches(patch_idx)%trace_types)
  deallocate(patches(patch_idx)%values)
end do
deallocate(patches)
deallocate(volumes(1)%window_shape)
deallocate(volumes)

end program example_9
