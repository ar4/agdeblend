#ifdef AGD_MPI
#include "mpi.h"
#endif /* AGD_MPI */
#include "agdeblend.h"

int deblend_c(int const n_patches, int const *const volumes,
              int const *const n_dims, int const *const *const window_shapes,
              int const *const *const coords, int const *const *const shapes,
              int long const *const *const shottimes,
              int const *const *const channels,
              enum AGDTraceType const *const *const trace_types,
              int const *const wavelet_lengths,
              int const *const *const wavelet_idxs,
              AGD_TYPE const *const *const wavelets,
              AGD_TYPE const initial_factor, int const n_its,
              int const print_freq,
#ifdef AGD_MPI
              MPI_Fint comm,
#endif /* AGD_MPI */
              AGD_TYPE *const *const values) {
  return agd_deblend(n_patches, volumes, n_dims, window_shapes, coords, shapes,
                     shottimes, channels, trace_types, wavelet_lengths,
                     wavelet_idxs, wavelets, initial_factor, n_its, print_freq,
#ifdef AGD_MPI
                     MPI_Comm_f2c(comm),
#endif /* AGD_MPI */
                     values);
}

int blend_c(int n_patches, int const *n_traces, int const *n_times,
            long int const *const *shottimes, int const *const *channels,
            enum AGDTraceType const *const *trace_types,
            AGD_TYPE *const *values, enum AGDBlendMode blend_mode,
            int taper_length, int n_patches_out, int const *n_traces_out,
            int const *n_times_out, long int const *const *shottimes_out,
            int const *const *channels_out,
            enum AGDTraceType const *const *trace_types_out,
#ifdef AGD_MPI
            MPI_Fint comm,
#endif /* AGD_MPI */
            AGD_TYPE *const *values_out) {
  return agd_blend(n_patches, n_traces, n_times,
                   shottimes, channels, trace_types, values,
                   blend_mode, taper_length, n_patches_out,
                   n_traces_out, n_times_out, shottimes_out,
                   channels_out, trace_types_out,
#ifdef AGD_MPI
                   MPI_Comm_f2c(comm),
#endif /* AGD_MPI */
                   values_out);
}
