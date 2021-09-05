#ifdef AGD_DOUBLE
#define C_AGD_TYPE c_double
#else
#define C_AGD_TYPE c_float
#endif /* AGD_DOUBLE */

module agdeblend_m
  use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float,        &
                                         c_double, c_ptr, c_loc
  use, intrinsic :: iso_fortran_env, only: error_unit
#ifdef AGD_MPI
#ifdef AGD_F08
  use mpi_f08
#else
  use mpi
#endif /* AGD_F08 */
#endif /* AGD_MPI*/
  implicit none

  private
  public :: blend, deblend, volume, deblend_patch, blend_patch,         &
            wavelet, AGDLive, AGDBad, AGDMissing, AGDBlendSum,          &
            AGDBlendMean, AGDBlendOverwrite

  enum, bind(c)
    enumerator :: AGDLive
    enumerator :: AGDBad
    enumerator :: AGDMissing
  end enum

  enum, bind(c)
    enumerator :: AGDBlendSum
    enumerator :: AGDBlendMean
    enumerator :: AGDBlendOverwrite
  end enum

  type volume
    integer(kind=c_int) :: n_dims
    integer(kind=c_int), allocatable :: window_shape(:)
  end type volume

  type deblend_patch
    integer(kind=c_int) :: volume_idx
    integer(kind=c_int), allocatable :: coords(:)
    integer(kind=c_int), allocatable :: data_shape(:)
    integer(kind=c_long), allocatable :: shottimes(:)
    integer(kind=c_int), allocatable :: channels(:)
    integer(kind=c_int), allocatable :: trace_types(:)
    integer(kind=c_int), allocatable :: wavelet_idxs(:)
    real(kind=C_AGD_TYPE), allocatable :: values(:)
  end type deblend_patch

  type blend_patch
    integer(kind=c_int) :: n_traces
    integer(kind=c_int) :: n_times
    integer(kind=c_long), allocatable :: shottimes(:)
    integer(kind=c_int), allocatable :: channels(:)
    integer(kind=c_int), allocatable :: trace_types(:)
    real(kind=C_AGD_TYPE), allocatable :: values(:)
  end type blend_patch

  type wavelet
    real(kind=C_AGD_TYPE), allocatable :: values(:)
  end type wavelet

  type row_major_array
    integer(kind=c_int), allocatable :: array(:)
  end type row_major_array

  interface

    function deblend_c(n_patches, volume_idxs, n_dims, window_shapes,   &
                       coords, shapes, shottimes, channels, trace_types,&
                       wavelet_lengths, wavelet_idxs, wavelets,         &
                       initial_factor, n_its, print_freq,               &
#ifdef AGD_MPI
                       comm,                                            &
#endif /* AGD_MPI */
                       values) bind(c)
      use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double,  &
                                             c_ptr
      implicit none
      integer(kind=c_int) :: deblend_c
      integer(kind=c_int), intent(in), value :: n_patches
      integer(kind=c_int), intent(in) :: volume_idxs(*)
      integer(kind=c_int), intent(in) :: n_dims(*)
      type(c_ptr), intent(in) :: window_shapes(*)
      type(c_ptr), intent(in) :: coords(*)
      type(c_ptr), intent(in) :: shapes(*)
      type(c_ptr), intent(in) :: shottimes(*)
      type(c_ptr), intent(in) :: channels(*)
      type(c_ptr), intent(in) :: trace_types(*)
      integer(kind=c_int), intent(in) :: wavelet_lengths(*)
      type(c_ptr), intent(in) :: wavelet_idxs(*)
      type(c_ptr), intent(in) :: wavelets(*)
      real(kind=C_AGD_TYPE), intent(in), value :: initial_factor
      integer(kind=c_int), intent(in), value :: n_its
      integer(kind=c_int), intent(in), value :: print_freq
#ifdef AGD_MPI
      integer, intent(in), value :: comm
#endif /* AGD_MPI */
      type(c_ptr), intent(in) :: values(*)
    end function deblend_c

    function blend_c(n_patches, n_traces, n_times, shottimes, channels, &
                     trace_types, values, blend_mode, taper_length,     &
                     n_patches_out, n_traces_out, n_times_out,          &
                     shottimes_out, channels_out, trace_types_out,      &
#ifdef AGD_MPI
                     comm,                                              &
#endif /* AGD_MPI */
                     values_out) bind(c)
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      integer(kind=c_int) :: blend_c
      integer(kind=c_int), intent(in), value :: n_patches
      integer(kind=c_int), intent(in) :: n_traces(*)
      integer(kind=c_int), intent(in) :: n_times(*)
      type(c_ptr), intent(in) :: shottimes(*)
      type(c_ptr), intent(in) :: channels(*)
      type(c_ptr), intent(in) :: trace_types(*)
      type(c_ptr), intent(in) :: values(*)
      integer(kind=c_int), intent(in), value :: blend_mode
      integer(kind=c_int), intent(in), value :: taper_length
      integer(kind=c_int), intent(in), value :: n_patches_out
      integer(kind=c_int), intent(in) :: n_traces_out(*)
      integer(kind=c_int), intent(in) :: n_times_out(*)
      type(c_ptr), intent(in) :: shottimes_out(*)
      type(c_ptr), intent(in) :: channels_out(*)
      type(c_ptr), intent(in) :: trace_types_out(*)
#ifdef AGD_MPI
      integer, intent(in), value :: comm
#endif /* AGD_MPI */
      type(c_ptr), intent(in) :: values_out(*)
    end function blend_c

  end interface

  contains

  subroutine deblend(patches, volumes, initial_factor, n_its,           &
                     print_freq,                                        &
#ifdef AGD_MPI
                     comm,                                              &
#endif /* AGD_MPI */
                     ierr, wavelets)

  type(deblend_patch), intent(in), allocatable, target :: patches(:)
  type(volume), intent(in), allocatable, target :: volumes(:)
  integer, intent(in) :: n_its
  real(kind=C_AGD_TYPE), intent(in) :: initial_factor
  integer, intent(in) :: print_freq
#ifdef AGD_MPI
#ifdef AGD_F08
  type(MPI_Comm), intent(in) :: comm
#else
  integer, intent(in) :: comm
#endif /* AGD_F08 */
#endif /* AGD_MPI */
  integer, intent(out) :: ierr
  type(wavelet), intent(in), allocatable, target, optional ::           &
    wavelets(:)

  integer(kind=c_int) :: n_patches
  integer(kind=c_int), allocatable :: volume_idxs(:)
  integer(kind=c_int), allocatable :: n_dims(:)
  type(row_major_array), allocatable, target :: window_shapes_row(:)
  type(row_major_array), allocatable, target :: coords_row(:)
  type(row_major_array), allocatable, target :: shapes_row(:)
  type(c_ptr), allocatable :: window_shapes(:)
  type(c_ptr), allocatable :: coords(:)
  type(c_ptr), allocatable :: shapes(:)
  type(c_ptr), allocatable :: shottimes(:)
  type(c_ptr), allocatable :: channels(:)
  type(c_ptr), allocatable :: trace_types(:)
  integer(kind=c_int), allocatable :: wavelet_lengths(:)
  type(c_ptr), allocatable :: wavelet_idxs(:)
  type(c_ptr), allocatable :: wavelet_values(:)
  type(c_ptr), allocatable :: values(:)
  integer :: array_n_dims
  integer :: n_volumes
  integer :: n_wavelets
  integer :: patch_idx
  integer :: volume_idx
  integer :: wavelet_idx
  integer :: idx

  ! Check inputs and get n_patches and n_wavelets
  if (.not. allocated(patches)) then
    write (error_unit, *) "ERROR: patches not allocated"
    ierr = 1
    return
  end if
  n_patches = size(patches)

  if (.not. allocated(volumes)) then
    write (error_unit, *) "ERROR: volumes not allocated"
    ierr = 1
    return
  end if
  n_volumes = size(volumes)

  if (.not. present(wavelets)) then
    n_wavelets = 0
  else if (.not. allocated(wavelets)) then
    n_wavelets = 0
  else
    n_wavelets = size(wavelets)
    do wavelet_idx = 1, n_wavelets
      if (.not. allocated(wavelets(wavelet_idx)%values)) then
        write (error_unit, *) "ERROR: wavelet values not allocated"
        ierr = 1
        return
      end if
    end do
  end if

  ! Convert from column-major to row-major
  allocate(window_shapes_row(n_volumes))
  allocate(coords_row(n_patches))
  allocate(shapes_row(n_patches))
  ! window_shapes
  do volume_idx = 1, n_volumes
    array_n_dims = volumes(volume_idx)%n_dims
    allocate(window_shapes_row(volume_idx)%array(array_n_dims))
    do idx = 1, array_n_dims
      window_shapes_row(volume_idx)%array(idx) =                        &
        volumes(volume_idx)%window_shape(array_n_dims - idx + 1)
    end do
  end do
  do patch_idx = 1, n_patches
    ! coords
    array_n_dims = volumes(patches(patch_idx)%volume_idx + 1)%n_dims
    allocate(coords_row(patch_idx)%array(array_n_dims - 1))
    do idx = 1, array_n_dims - 1
      coords_row(patch_idx)%array(idx) =                                &
        patches(patch_idx)%coords(array_n_dims - 1 - idx + 1)
    end do
    ! data_shapes
    allocate(shapes_row(patch_idx)%array(array_n_dims))
    do idx = 1, array_n_dims
      shapes_row(patch_idx)%array(idx) =                                &
        patches(patch_idx)%data_shape(array_n_dims - idx + 1)
    end do
  end do

  ! The C code requires a different data structure, so we need
  ! to allocate and set new arrays
  allocate(volume_idxs(n_patches))
  allocate(n_dims(n_volumes))
  allocate(window_shapes(n_volumes))
  allocate(coords(n_patches))
  allocate(shapes(n_patches))
  allocate(shottimes(n_patches))
  allocate(channels(n_patches))
  allocate(trace_types(n_patches))
  if (n_wavelets > 0) then
    allocate(wavelet_lengths(n_wavelets))
    allocate(wavelet_idxs(n_patches))
    allocate(wavelet_values(n_wavelets))
  end if
  allocate(values(n_patches))

  ! In most cases, we only need to pass an array of pointers,
  ! so we do not need to copy the data
  do patch_idx = 1, n_patches
    volume_idxs(patch_idx) = patches(patch_idx)%volume_idx
    coords(patch_idx) = c_loc(coords_row(patch_idx)%array)
    shapes(patch_idx) = c_loc(shapes_row(patch_idx)%array)
    shottimes(patch_idx) = c_loc(patches(patch_idx)%shottimes)
    channels(patch_idx) = c_loc(patches(patch_idx)%channels)
    trace_types(patch_idx) = c_loc(patches(patch_idx)%trace_types)
    if (n_wavelets > 0) then
      wavelet_idxs(patch_idx) = c_loc(patches(patch_idx)%wavelet_idxs)
    end if
    values(patch_idx) = c_loc(patches(patch_idx)%values)
  end do

  do volume_idx = 1, n_volumes
    n_dims(volume_idx) = volumes(volume_idx)%n_dims
    window_shapes(volume_idx) =                                         &
      c_loc(window_shapes_row(volume_idx)%array)
  end do

  if (n_wavelets > 0) then
    do wavelet_idx = 1, n_wavelets
      wavelet_lengths(wavelet_idx) = size(wavelets(wavelet_idx)%values)
      wavelet_values(wavelet_idx) = c_loc(wavelets(wavelet_idx)%values)
    end do
  end if

  ierr = deblend_c(n_patches, volume_idxs, n_dims, window_shapes,       &
                   coords, shapes, shottimes, channels, trace_types,    &
                   wavelet_lengths, wavelet_idxs, wavelet_values,       &
                   initial_factor, n_its, print_freq,                   &
#ifdef AGD_MPI
#ifdef AGD_F08
                   comm%mpi_val,                                        &
#else
                   comm,                                                &
#endif /* AGD_F08 */
#endif /* AGD_MPI */
                   values)

  ! Deallocate row-major arrays
  ! window_shapes
  do volume_idx = 1, n_volumes
    deallocate(window_shapes_row(volume_idx)%array)
  end do
  do patch_idx = 1, n_patches
    ! coords
    deallocate(coords_row(patch_idx)%array)
    ! data_shapes
    deallocate(shapes_row(patch_idx)%array)
  end do
  deallocate(window_shapes_row)
  deallocate(coords_row)
  deallocate(shapes_row)

  deallocate(volume_idxs)
  deallocate(n_dims)
  deallocate(window_shapes)
  deallocate(coords)
  deallocate(shapes)
  deallocate(shottimes)
  deallocate(channels)
  deallocate(trace_types)
  if (n_wavelets > 0) then
    deallocate(wavelet_lengths)
    deallocate(wavelet_idxs)
    deallocate(wavelet_values)
  end if
  deallocate(values)

  end subroutine deblend

  subroutine blend(patches, blend_mode,                                 &
#ifdef AGD_MPI
                   comm,                                                &
#endif /* AGD_MPI */
                   ierr, taper_length, patches_out)

  type(blend_patch), intent(in), allocatable, target :: patches(:)
  integer(kind=c_int), intent(in) :: blend_mode
#ifdef AGD_MPI
#ifdef AGD_F08
  type(MPI_Comm), intent(in) :: comm
#else
  integer, intent(in) :: comm
#endif /* AGD_F08 */
#endif /* AGD_MPI */
  integer, intent(out) :: ierr
  integer(kind=c_int), intent(in), optional :: taper_length
  type(blend_patch), intent(in), allocatable, target, optional ::       &
    patches_out(:)

  integer(kind=c_int) :: n_patches
  integer(kind=c_int) :: n_patches_out
  integer(kind=c_int), allocatable :: n_traces(:)
  integer(kind=c_int), allocatable :: n_times(:)
  type(c_ptr), allocatable :: shottimes(:)
  type(c_ptr), allocatable :: channels(:)
  type(c_ptr), allocatable :: trace_types(:)
  type(c_ptr), allocatable :: values(:)
  integer(kind=c_int), allocatable :: n_traces_out(:)
  integer(kind=c_int), allocatable :: n_times_out(:)
  type(c_ptr), allocatable :: shottimes_out(:)
  type(c_ptr), allocatable :: channels_out(:)
  type(c_ptr), allocatable :: trace_types_out(:)
  type(c_ptr), allocatable :: values_out(:)
  integer(kind=c_int) :: taper_length_o
  integer :: patch_idx

  if (.not. allocated(patches)) then
    write (error_unit, *) "ERROR: patches not allocated"
    ierr = 1
    return
  end if
  n_patches = size(patches)

  if (.not. present(patches_out)) then
    n_patches_out = n_patches
  else if (.not. allocated(patches_out)) then
    n_patches_out = n_patches
  else
    n_patches_out = size(patches_out)
  end if

  if (blend_mode == AGDBlendMean .and. .not. present(taper_length)) then
    write (error_unit, *) "ERROR: taper_length is required with " //    &
                          "AGDBlendMean"
    ierr = 1
    return
  end if

  if (present(taper_length)) then
    taper_length_o = taper_length
  else
    taper_length_o = 0
  end if

  allocate(n_traces(n_patches))
  allocate(n_times(n_patches))
  allocate(shottimes(n_patches))
  allocate(channels(n_patches))
  allocate(trace_types(n_patches))
  allocate(values(n_patches))
  allocate(n_traces_out(n_patches_out))
  allocate(n_times_out(n_patches_out))
  allocate(shottimes_out(n_patches_out))
  allocate(channels_out(n_patches_out))
  allocate(trace_types_out(n_patches_out))
  allocate(values_out(n_patches_out))

  do patch_idx = 1, n_patches
    n_traces(patch_idx) = patches(patch_idx)%n_traces
    n_times(patch_idx) = patches(patch_idx)%n_times
    shottimes(patch_idx) = c_loc(patches(patch_idx)%shottimes)
    channels(patch_idx) = c_loc(patches(patch_idx)%channels)
    trace_types(patch_idx) = c_loc(patches(patch_idx)%trace_types)
    values(patch_idx) = c_loc(patches(patch_idx)%values)
  end do

  do patch_idx = 1, n_patches_out
    if (.not. present(patches_out)) then
      ! patches_out not provided, so reuse patches
      n_traces_out(patch_idx) = patches(patch_idx)%n_traces
      n_times_out(patch_idx) = patches(patch_idx)%n_times
      shottimes_out(patch_idx) = c_loc(patches(patch_idx)%shottimes)
      channels_out(patch_idx) = c_loc(patches(patch_idx)%channels)
      trace_types_out(patch_idx) = c_loc(patches(patch_idx)%trace_types)
      values_out(patch_idx) = c_loc(patches(patch_idx)%values)
    else
      n_traces_out(patch_idx) = patches_out(patch_idx)%n_traces
      n_times_out(patch_idx) = patches_out(patch_idx)%n_times
      shottimes_out(patch_idx) = c_loc(patches_out(patch_idx)%shottimes)
      channels_out(patch_idx) = c_loc(patches_out(patch_idx)%channels)
      trace_types_out(patch_idx) =                                      &
        c_loc(patches_out(patch_idx)%trace_types)
      if (allocated(patches_out(patch_idx)%values)) then
        values_out(patch_idx) = c_loc(patches_out(patch_idx)%values)
      else
        values_out(patch_idx) = c_loc(patches(patch_idx)%values)
      end if
    end if
  end do

  ierr = blend_c(n_patches, n_traces, n_times, shottimes, channels,     &
                 trace_types, values, blend_mode, taper_length_o,       &
                 n_patches_out, n_traces_out, n_times_out,              &
                 shottimes_out, channels_out, trace_types_out,          &
#ifdef AGD_MPI
#ifdef AGD_F08
                 comm%mpi_val,                                          &
#else
                 comm,                                                  &
#endif /* AGD_F08 */
#endif /* AGD_MPI */
                 values_out)

  deallocate(n_traces)
  deallocate(n_times)
  deallocate(shottimes)
  deallocate(channels)
  deallocate(trace_types)
  deallocate(values)
  deallocate(n_traces_out)
  deallocate(n_times_out)
  deallocate(shottimes_out)
  deallocate(channels_out)
  deallocate(trace_types_out)
  deallocate(values_out)

  end subroutine blend

end module agdeblend_m
