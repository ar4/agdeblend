#ifndef AGD_H
#define AGD_H

#ifdef AGD_MPI
#include <mpi.h>
#endif /* AGD_MPI */

enum AGDTraceType { AGDLive, AGDBad, AGDMissing };
enum AGDBlendMode { AGDBlendSum, AGDBlendMean, AGDBlendOverwrite };

#ifdef AGD_DOUBLE
#define AGD_TYPE double
#else
#define AGD_TYPE float
#endif

#define AGD_PAD_FRACTION ((AGD_TYPE)0.25)
#define AGD_TARGET_LAMB ((AGD_TYPE)0.001) /* 1e-3 */
#define AGD_MAX_N_COORDS 4
#define AGD_FFTW_FLAGS FFTW_ESTIMATE

int agd_deblend(int n_patches, int const *volumes, int const *n_dims,
                int const *const *window_shapes, int const *const *coords,
                int const *const *shapes, long int const *const *shottimes,
                int const *const *channels,
                enum AGDTraceType const *const *trace_types,
                int const *wavelet_lengths, int const *const *wavelet_idxs,
                AGD_TYPE const *const *wavelets, AGD_TYPE const initial_factor, int n_its, int print_freq,
#ifdef AGD_MPI
                MPI_Comm comm,
#endif /* AGD_MPI */
                AGD_TYPE *const *data);
int agd_blend(int n_patches, int const *n_traces, int const *n_times,
              long int const *const *shottimes, int const *const *channels,
              enum AGDTraceType const *const *trace_types,
              AGD_TYPE *const *data, enum AGDBlendMode blend_mode,
              int taper_length, int n_patches_out, int const *n_traces_out,
              int const *n_times_out, long int const *const *shottimes_out,
              int const *const *channels_out,
              enum AGDTraceType const *const *trace_types_out,
#ifdef AGD_MPI
              MPI_Comm comm,
#endif /* AGD_MPI */
              AGD_TYPE *const *data_out);

#endif
