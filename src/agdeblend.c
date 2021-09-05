/*
AGDeblend
Copyright (C) 2021 Alan Richardson

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef AGD_MPI
#include <mpi.h>
#endif /* AGD_MPI */
#ifdef AGD_THREADS
#include <omp.h>
#endif /* AGD_THREADS */

#include "agdeblend.h"

/* Defines */

#ifdef __cplusplus
#define ZERO_INIT \
  {}
#else
#define ZERO_INIT \
  { 0 }
#endif
#ifndef M_PI
#define M_PI 3.14159265359
#endif
#ifdef AGD_DOUBLE
#define AGD_MPI_TYPE MPI_DOUBLE
#define AGD_ZERO 0.0
#define AGD_ONE 1.0
#define AGD_TWO 2.0
#define AGD_HALF 0.5
#define AGD_COS cos
#define AGD_SIN sin
#define AGD_ATAN2 atan2
#define AGD_EXP exp
#define AGD_LOG log
#define AGD_POW pow
#define AGD_SQRT sqrt
#define AGD_CEIL ceil
#define AGD_FFTW_PLAN fftw_plan
#define AGD_FFTW_COMPLEX fftw_complex
#define AGD_FFTW_PLAN_C2R fftw_plan_many_dft_c2r
#define AGD_FFTW_PLAN_R2C fftw_plan_many_dft_r2c
#define AGD_FFTW_DESTROY_PLAN fftw_destroy_plan
#define AGD_FFTW_EXECUTE fftw_execute
#define AGD_FFTW_EXECUTE_C2R fftw_execute_dft_c2r
#define AGD_FFTW_EXECUTE_R2C fftw_execute_dft_r2c
#define AGD_FFTW_ALLOC_COMPLEX fftw_alloc_complex
#define AGD_FFTW_MALLOC fftw_malloc
#define AGD_FFTW_FREE fftw_free
#else
#define AGD_MPI_TYPE MPI_FLOAT
#define AGD_ZERO 0.0f
#define AGD_ONE 1.0f
#define AGD_TWO 2.0f
#define AGD_HALF 0.5f
#define AGD_COS cosf
#define AGD_SIN sinf
#define AGD_ATAN2 atan2f
#define AGD_EXP expf
#define AGD_LOG logf
#define AGD_POW powf
#define AGD_SQRT sqrtf
#define AGD_CEIL ceilf
#define AGD_FFTW_PLAN fftwf_plan
#define AGD_FFTW_COMPLEX fftwf_complex
#define AGD_FFTW_PLAN_C2R fftwf_plan_many_dft_c2r
#define AGD_FFTW_PLAN_R2C fftwf_plan_many_dft_r2c
#define AGD_FFTW_DESTROY_PLAN fftwf_destroy_plan
#define AGD_FFTW_EXECUTE fftwf_execute
#define AGD_FFTW_EXECUTE_C2R fftwf_execute_dft_c2r
#define AGD_FFTW_EXECUTE_R2C fftwf_execute_dft_r2c
#define AGD_FFTW_ALLOC_COMPLEX fftwf_alloc_complex
#define AGD_FFTW_MALLOC fftwf_malloc
#define AGD_FFTW_FREE fftwf_free
#endif

/* Data structures */

struct SpaceWindow {
  int *coords;
  int *space_shape;
  int *taper_length;
  int *use_taper;
  int patch_idx;
  int n_dims;
  int n_traces;
};

struct WindowArrays {
  long int *shottimes;
  int *channels;
  enum AGDTraceType *trace_types;
  int *wavelet_idxs;
};

struct FKWindowConfig {
  AGD_FFTW_PLAN fwd_plan;
  AGD_FFTW_PLAN adj_plan;
  int *window_shape;
  int *padded_window_shape;
  int *padded_window_shape_fk;
  int *use_taper;
  AGD_TYPE **taper_windows;
  AGD_TYPE scaler;
  int n_time_windows;
  int n_dims;
  int n_times;
  int n_times_output_capacity;
  int n_traces;
  int n_traces_padded;
  int n_elements_tx;
  int n_elements_tx_output_capacity;
  int n_elements_padded_window;
  int n_elements_padded_window_fk;
};

struct FKConfigs {
  struct FKWindowConfig *window_configs;
  int n_window_configs;
};

struct WaveletParams {
  int *source_config_idxs;
  int plan_idx;
};

struct WaveletSourceConfig {
  AGD_FFTW_COMPLEX *wavelet;
  int original_idx;
  int n_times;
};

struct WaveletFFTWPlan {
  AGD_FFTW_PLAN forward_plan;
  AGD_FFTW_PLAN backward_plan;
  int n_traces;
  int n_times;
};

struct WaveletConfig {
  struct WaveletSourceConfig *source_configs;
  struct WaveletFFTWPlan *plans;
  int n_source_configs;
  int n_fftw_plans;
};

struct Interval {
  long int start;
  long int stop;
};

struct BlendCoords {
  size_t interval_start;
  int trace_start;
  int channel_idx;
  int interval_idx;
  int n_times;
};

struct ChannelIntervals {
  struct Interval *intervals;
  int n_intervals;
};

struct BlendParams {
  struct BlendCoords *blend_coords;
  int n_traces;
  int trace_length;
};

struct BlendConfig {
  int *unique_channel_values;
  struct ChannelIntervals *channels_intervals;
  struct BlendCoords *mute_coords;
#ifdef AGD_MPI
  struct BlendCoords **ranks_coords;
  AGD_TYPE **send_buffers;
  AGD_TYPE **receive_buffers;
  MPI_Request *requests;
  MPI_Comm comm;
  int *n_overlaps;
  int comm_size;
#endif /* AGD_MPI */
  int n_channels;
  int n_mute_coords;
};

struct ModelWindow {
  struct SpaceWindow space_window;
  struct WaveletParams wavelet_params;
  struct BlendParams blend_params;
  int fk_config_idx;
};

struct Model {
  struct FKConfigs fk_configs;
  struct WaveletConfig wavelet_config;
  struct BlendConfig blend_config;
  struct ModelWindow *windows;
  int *n_times_fk;
  int *n_times_conv;
  void **temp_1;
  void **temp_2;
  int n_windows;
};

struct PatchProperties {
  int *n_traces;
  int *n_times_data;
  int *n_times_window;
  int *n_times_wavelet;
};

struct VolumeCoords {
  int patch_idx;
  int *coords;
  int *shape;
};

#include "utils.c"
#include "fk.c"
#include "blend.c"
#include "wavelet.c"
#include "model.c"

static int check_positive(int const value, char const name[]) {
  if (value <= 0) {
    fprintf(stderr, "ERROR: %s must be > 0\n", name);
    return 1;
  }
  return 0;
}

static int check_non_null(void const *const ptr, char const name[]) {
  if (ptr == NULL) {
    fprintf(stderr, "ERROR: %s must not be NULL\n", name);
    return 1;
  }
  return 0;
}

static int check_deblend_inputs(
    int const n_patches, int const *const volumes, int const *const n_dims,
    int const *const *const window_shapes, int const *const *const coords,
    int const *const *const shapes, long int const *const *const shottimes,
    int const *const *const channels,
    enum AGDTraceType const *const *const trace_types,
    int const *const wavelet_lengths, int const *const *const wavelet_idxs,
    AGD_TYPE const *const *const wavelets, AGD_TYPE const initial_factor,
    int const n_its,
#ifdef AGD_MPI
    MPI_Comm comm,
#endif /* AGD_MPI */
    AGD_TYPE const *const *const data) {
  int patch_idx;
#ifdef AGD_MPI
  int comm_rank;
  int comm_size;
  int rank_idx;
  int n_volumes = get_n_volumes(n_patches, volumes);
#endif

  if (check_positive(n_patches, "n_patches")) goto err;
  if (check_positive(n_its, "n_its")) goto err;

  if (initial_factor <= AGD_ZERO || initial_factor > AGD_ONE) {
    fprintf(stderr, "ERROR: initial_factor must be in (0, 1]\n");
    goto err;
  }

  if (check_non_null(volumes, "volumes")) goto err;
  if (check_non_null(n_dims, "n_dims")) goto err;
  if (check_non_null(window_shapes, "window_shapes")) goto err;
  if (check_non_null(coords, "coords")) goto err;
  if (check_non_null(shapes, "shapes")) goto err;
  if (check_non_null(shottimes, "shottimes")) goto err;
  if (check_non_null(channels, "channels")) goto err;
  if (check_non_null(trace_types, "trace_types")) goto err;
  if (check_non_null(data, "data")) goto err;

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const volume_idx = volumes[patch_idx];
    size_t window_size = 1;
    size_t n_traces = 1;
    int max_n_times = shapes[patch_idx][n_dims[volume_idx] - 1] +
                      window_shapes[volume_idx][n_dims[volume_idx] - 1];
    int dim_idx;
    if (check_non_null(coords[patch_idx], "coords")) goto err;
    if (check_non_null(shottimes[patch_idx], "shottimes")) goto err;
    if (check_non_null(channels[patch_idx], "channels")) goto err;
    if (check_non_null(trace_types[patch_idx], "trace_types")) goto err;
    if (check_non_null(data[patch_idx], "data")) goto err;

    if (n_dims[volume_idx] <= 1) {
      fprintf(stderr, "ERROR: n_dims must be > 1\n");
      goto err;
    }
    if (n_dims[volume_idx] - 1 > AGD_MAX_N_COORDS) {
      fprintf(
          stderr,
          "ERROR: AGD_MAX_N_COORDS is currently set to %d, "
          "but %d space dimensions have been requested. "
          "Please increase AGD_MAX_N_COORDS in the source code and recompile\n",
          AGD_MAX_N_COORDS, n_dims[volume_idx] - 1);
      goto err;
    }
    if (check_non_null(window_shapes[volume_idx], "window_shapes")) goto err;
    if (check_non_null(shapes[patch_idx], "shapes")) goto err;
    for (dim_idx = 0; dim_idx < n_dims[volume_idx]; dim_idx++) {
      if (window_shapes[volume_idx][dim_idx] < 2) {
        fprintf(stderr, "ERROR: window_shapes must be at least 2\n");
        goto err;
      }
      if ((window_shapes[volume_idx][dim_idx] != shapes[patch_idx][dim_idx]) &&
          (window_shapes[volume_idx][dim_idx] % 2 != 0)) {
        fprintf(stderr,
                "ERROR: window_shapes must be even if it doesn't "
                "cover the whole dimension\n");
        goto err;
      }
      if (window_shapes[volume_idx][dim_idx] > shapes[patch_idx][dim_idx]) {
        fprintf(stderr,
                "ERROR: window_shape must not be bigger than patch shape\n");
        goto err;
      }
    }
    for (dim_idx = 0; dim_idx < n_dims[volume_idx] - 1; dim_idx++) {
      window_size *= (size_t)window_shapes[volume_idx][dim_idx];
      n_traces *= (size_t)shapes[patch_idx][dim_idx];
    }
    if (2 * window_size * (size_t)max_n_times >= INT_MAX) {
      fprintf(stderr, "ERROR: The window for volume %d is too big\n",
              volume_idx);
      goto err;
    }
    if (n_traces >= INT_MAX) {
      fprintf(stderr, "ERROR: Patches must have fewer than %d traces\n",
              INT_MAX);
      goto err;
    }
  }

  if (wavelet_lengths != NULL || wavelet_idxs != NULL || wavelets != NULL) {
    int wavelet_idx;
    int n_wavelets = 0;
    int *n_traces = NULL;
    int *volume_wavelet_lengths = NULL;
    if (check_non_null(wavelet_lengths, "wavelet_lengths")) goto err;
    if (check_non_null(wavelet_idxs, "wavelet_idxs")) goto err;
    if (check_non_null(wavelets, "wavelets")) goto err;

    if (set_n_traces_per_patch(n_patches, volumes, n_dims, shapes, &n_traces))
      goto err;

    if (set_volume_wavelet_lengths(n_patches, volumes, n_traces, trace_types,
                                   wavelet_lengths, wavelet_idxs,
                                   &volume_wavelet_lengths)) {
      free(volume_wavelet_lengths);
      free(n_traces);
      goto err;
    }

    for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      int const volume_idx = volumes[patch_idx];
      if (check_non_null(wavelet_idxs[patch_idx], "wavelet_idxs")) {
        free(volume_wavelet_lengths);
        free(n_traces);
        goto err;
      }
      if (n_traces[patch_idx] > 0) {
        int trace_idx;
        if (volume_wavelet_lengths[volume_idx] == 0) {
          /* Check that patch contains at least one live trace */
          for (trace_idx = 0; trace_idx < n_traces[patch_idx]; trace_idx++) {
            if (trace_types[patch_idx][trace_idx] == AGDLive) {
              break;
            }
          }
          if (trace_idx == n_traces[patch_idx])
            continue; /* No live traces in patch */
        }
        for (trace_idx = 0; trace_idx < n_traces[patch_idx]; trace_idx++) {
          if (trace_types[patch_idx][trace_idx] != AGDLive) continue;
          wavelet_idx = wavelet_idxs[patch_idx][trace_idx];
          if ((wavelet_idx < 0)) {
            fprintf(stderr, "ERROR: wavelet_idx must not be negative\n");
            free(volume_wavelet_lengths);
            free(n_traces);
            goto err;
          }
          if (wavelet_idx + 1 > n_wavelets) n_wavelets = wavelet_idx + 1;
          if (wavelet_lengths[wavelet_idxs[patch_idx][trace_idx]] !=
              volume_wavelet_lengths[volume_idx]) {
            fprintf(stderr,
                    "ERROR: wavelets in the same volume must have the same "
                    "length\n");
            free(volume_wavelet_lengths);
            free(n_traces);
            goto err;
          }
        }
      }
    }
    free(volume_wavelet_lengths);
    free(n_traces);

    for (wavelet_idx = 0; wavelet_idx < n_wavelets; wavelet_idx++) {
      if (check_non_null(wavelets[wavelet_idx], "wavelets")) goto err;
    }
  }

#ifdef AGD_MPI
  if (MPI_Comm_rank(comm, &comm_rank)) goto err;
  if (MPI_Comm_size(comm, &comm_size)) goto err;
  if (MPI_Allreduce(MPI_IN_PLACE, &n_volumes, 1, MPI_INT, MPI_MAX, comm))
    goto err;

  for (rank_idx = 1; rank_idx < comm_size; rank_idx++) {
    if (comm_rank == 0) {
      int rank_n_its;
      int *rank_n_dims = (int *)malloc((size_t)n_volumes * sizeof(int));
      int volume_idx;

      if (rank_n_dims == NULL) goto err;

      if (MPI_Recv(&rank_n_its, 1, MPI_INT, rank_idx, 0, comm,
                   MPI_STATUS_IGNORE))
        goto err;
      if (rank_n_its != n_its) {
        fprintf(stderr, "ERROR: n_its must be the same on all ranks\n");
        free(rank_n_dims);
        goto err;
      }

      if (MPI_Recv(rank_n_dims, n_volumes, MPI_INT, rank_idx, 0, comm,
                   MPI_STATUS_IGNORE))
        goto err;
      for (volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
        if (rank_n_dims[volume_idx] != n_dims[volume_idx]) {
          fprintf(stderr, "ERROR: n_dims must be the same on all ranks\n");
          free(rank_n_dims);
          goto err;
        }
      }

      for (volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
        int *rank_window_shape =
            (int *)malloc((size_t)n_dims[volume_idx] * sizeof(int));
        int dim_idx;
        if (rank_window_shape == NULL) goto err;
        if (MPI_Recv(rank_window_shape, n_dims[volume_idx], MPI_INT, rank_idx,
                     0, comm, MPI_STATUS_IGNORE))
          goto err;
        for (dim_idx = 0; dim_idx < n_dims[volume_idx]; dim_idx++) {
          if (rank_window_shape[dim_idx] !=
              window_shapes[volume_idx][dim_idx]) {
            fprintf(stderr,
                    "ERROR: window_shapes must be the same on all ranks\n");
            free(rank_n_dims);
            free(rank_window_shape);
            goto err;
          }
        }
        free(rank_window_shape);
        rank_window_shape = NULL;
      }

      free(rank_n_dims);
      rank_n_dims = NULL;

    } else if (comm_rank == rank_idx) {
      int volume_idx;
      if (MPI_Send(&n_its, 1, MPI_INT, 0, 0, comm)) goto err;
      if (MPI_Send(n_dims, n_volumes, MPI_INT, 0, 0, comm)) goto err;
      for (volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
        if (MPI_Send(window_shapes[volume_idx], n_dims[volume_idx], MPI_INT, 0,
                     0, comm))
          goto err;
      }
    }
  }
#endif /* AGD_MPI */

  return 0;
err:
  fprintf(stderr, "ERROR in agd_deblend inputs\n");
  return 1;
}

int agd_deblend(int const n_patches, int const *const volumes,
                int const *const n_dims, int const *const *const window_shapes,
                int const *const *const coords, int const *const *const shapes,
                long int const *const *const shottimes,
                int const *const *const channels,
                enum AGDTraceType const *const *const trace_types,
                int const *const wavelet_lengths,
                int const *const *const wavelet_idxs,
                AGD_TYPE const *const *const wavelets,
                AGD_TYPE const initial_factor, int const n_its,
                int const print_freq,
#ifdef AGD_MPI
                MPI_Comm comm,
#endif /* AGD_MPI */
                AGD_TYPE *const *const data) {
  AGD_TYPE const decay =
      AGD_EXP(AGD_LOG(AGD_TARGET_LAMB / initial_factor) / (AGD_TYPE)n_its);
  int const first_linear_it =
      (int)((AGD_LOG(decay) * ((AGD_TYPE)n_its - AGD_ONE) + AGD_ONE) /
            AGD_LOG(decay));
  AGD_TYPE linear_slope = initial_factor *
                          AGD_POW(decay, (AGD_TYPE)first_linear_it) /
                          (AGD_TYPE)(n_its - 1 - first_linear_it);
  AGD_TYPE stepsize;
  AGD_TYPE lamb;
  AGD_TYPE ***model_blended = NULL;
  AGD_FFTW_COMPLEX **model_fk = NULL;
  AGD_TYPE ***blended = NULL;
  struct Model model = ZERO_INIT;
  int it;
  int comm_rank = 0;

#ifdef AGD_MPI
  MPI_Comm_rank(comm, &comm_rank);
#endif /* AGD_MPI */

  /* Check the inputs */
  if (check_deblend_inputs(n_patches, volumes, n_dims, window_shapes, coords,
                           shapes, shottimes, channels, trace_types,
                           wavelet_lengths, wavelet_idxs, wavelets,
                           initial_factor, n_its,
#ifdef AGD_MPI
                           comm,
#endif /* AGD_MPI */
                           (AGD_TYPE const *const *)data))
    goto err;

  /* Get the model configuration */
  if (set_model(n_patches, volumes, n_dims, window_shapes, coords, shapes,
                shottimes, channels, trace_types, wavelet_lengths, wavelet_idxs,
                wavelets,
#ifdef AGD_MPI
                comm,
#endif /* AGD_MPI */
                &model))
    goto err;

  /* Allocate memory for the model in the FK and blended domains */
  if (allocate_model_fk(&model, &model_fk)) goto err;
  if (allocate_blended(&(model.blend_config), &model_blended)) goto err;

  /* Get the stepsize */
  stepsize = get_stepsize(&model, (AGD_TYPE *const *const *)model_blended);
  if (stepsize > (AGD_TYPE)0.3) {
    printf("WARNING: There does not appear to be any overlap between shots\n");
  }

  /* Blend the input data */
  if (blend_input((AGD_TYPE const *const *)data, n_patches, volumes, n_dims,
                  shapes, shottimes, channels, trace_types,
                  &(model.blend_config), &blended))
    goto err;

  /* Get initial model update and lamb value */
  model_adjoint_update((AGD_TYPE const *const *const *)blended, &model,
                       stepsize, model_fk);
  lamb = get_initial_lamb(&model, (AGD_FFTW_COMPLEX const *const *)model_fk);
  linear_slope *= lamb;
  lamb *= initial_factor;
  shrink_threshold(&model, lamb, model_fk);

  if (print_freq > 0) {
    AGD_TYPE const residual = get_residual_norm(
        &(model.blend_config), (AGD_TYPE const *const *const *)blended);
    if (comm_rank == 0) printf("Iteration 0 Residual %f\n", (double)residual);
  }

  /* Iteratively deblend using ISTA with a decaying threshold */
  for (it = 1; it < n_its; it++) {
    int print_residual = print_freq > 0 && it % print_freq == 0;

    /* y = A(x) */
    if (model_forward((AGD_FFTW_COMPLEX const *const *)model_fk, &model,
                      (AGD_TYPE *const *const *)model_blended))
      goto err;

    /* r = y - d */
    set_residual((AGD_TYPE const *const *const *)blended, &(model.blend_config),
                 (AGD_TYPE *const *const *)model_blended);

    if (print_residual) {
      AGD_TYPE const residual = get_residual_norm(
          &(model.blend_config), (AGD_TYPE const *const *const *)model_blended);
      if (comm_rank == 0)
        printf("Iteration %d Residual %f\n", it, (double)residual);
    }

    /* x = x - stepsize * A^H(r) */
    model_adjoint_update((AGD_TYPE const *const *const *)model_blended, &model,
                         stepsize, model_fk);

    /* Soft threshold x, shrinking amplitude by up to lamb */
    if (it < first_linear_it) {
      lamb *= decay;
    } else {
      lamb = (AGD_TYPE)(n_its - 1 - it) * linear_slope;
    }
    shrink_threshold(&model, lamb, model_fk);
  }

  /* Get the model output in the TX domain */
  set_output((AGD_FFTW_COMPLEX const *const *)model_fk, n_patches, volumes,
             n_dims, shapes, &model, data);

  /* Free the allocated memory and return */
  free_blended(&(model.blend_config), &blended);
  free_blended(&(model.blend_config), &model_blended);
  free_model_fk(model.n_windows, &model_fk);
  free_model(&model);
  return 0;

err:
  fprintf(stderr, "ERROR in agd_deblend\n");
  free_blended(&(model.blend_config), &blended);
  free_blended(&(model.blend_config), &model_blended);
  free_model_fk(model.n_windows, &model_fk);
  free_model(&model);
  return 1;
}

static int check_blend_inputs(
    int const n_patches, int const *const n_traces_per_volume,
    int const *const n_times_per_volume, long int const *const *const shottimes,
    int const *const *const channels,
    enum AGDTraceType const *const *const trace_types,
    AGD_TYPE *const *const data, enum AGDBlendMode const blend_mode,
    int const taper_length, int const n_patches_out,
    int const *const n_traces_per_volume_out,
    int const *const n_times_per_volume_out,
    long int const *const *const shottimes_out,
    int const *const *const channels_out,
    enum AGDTraceType const *const *const trace_types_out,
    AGD_TYPE *const *const data_out) {
  int patch_idx;

  if (check_positive(n_patches, "n_patches")) goto err;
  if (check_positive(n_patches, "n_patches_out")) goto err;

  if (check_non_null(n_traces_per_volume, "n_traces_per_volume")) goto err;
  if (check_non_null(n_traces_per_volume_out, "n_traces_per_volume_out"))
    goto err;
  if (check_non_null(n_times_per_volume, "n_times_per_volume")) goto err;
  if (check_non_null(n_times_per_volume_out, "n_times_per_volume_out"))
    goto err;
  if (check_non_null(shottimes, "shottimes")) goto err;
  if (check_non_null(shottimes_out, "shottimes_out")) goto err;
  if (check_non_null(channels, "channels")) goto err;
  if (check_non_null(channels_out, "channels_out")) goto err;
  if (check_non_null(trace_types, "trace_types")) goto err;
  if (check_non_null(trace_types_out, "trace_types_out")) goto err;
  if (check_non_null(data, "data")) goto err;
  if (check_non_null(data_out, "data_out")) goto err;

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    if (check_positive(n_traces_per_volume[patch_idx], "n_traces_per_volume"))
      goto err;
    if (check_positive(n_times_per_volume[patch_idx], "n_times_per_volume"))
      goto err;

    if (check_non_null(shottimes[patch_idx], "shottimes")) goto err;
    if (check_non_null(channels[patch_idx], "channels")) goto err;
    if (check_non_null(trace_types[patch_idx], "trace_types")) goto err;
    if (check_non_null(data[patch_idx], "data")) goto err;
  }

  for (patch_idx = 0; patch_idx < n_patches_out; patch_idx++) {
    int trace_idx;
    if (check_positive(n_traces_per_volume_out[patch_idx],
                       "n_traces_per_volume_out"))
      goto err;
    if (check_positive(n_times_per_volume_out[patch_idx],
                       "n_times_per_volume_out"))
      goto err;

    if (check_non_null(shottimes_out[patch_idx], "shottimes_out")) goto err;
    if (check_non_null(channels_out[patch_idx], "channels_out")) goto err;
    if (check_non_null(trace_types_out[patch_idx], "trace_types_out")) goto err;
    if (check_non_null(data_out[patch_idx], "data_out")) goto err;

    for (trace_idx = 0; trace_idx < n_traces_per_volume_out[patch_idx];
         trace_idx++) {
      if (trace_types_out[patch_idx][trace_idx] == AGDBad) {
        fprintf(stderr, "ERROR: Output trace types may not be AGDBad\n");
        goto err;
      }
    }
  }

  if (blend_mode == AGDBlendMean) {
    if (taper_length < 0) {
      fprintf(
          stderr,
          "ERROR: taper_length must be >= 0 when blend_mode is AGDBlendMean\n");
      goto err;
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in agd_blend inputs\n");
  return 1;
}

#ifdef AGD_MPI
/* When using MPI, we need to create blend_config using patches that cover
 * both the input and output. This function joins the input and output
 * patches so that this can be done.
 *
 * An example demonstrating why this is necessary is if on this rank
 * the input only covers times 0-8s, and on another rank the input
 * covers times 8-16s, but on this rank the requested output
 * covers times 8-16s. If only the input data was used when
 * creating blend_config, then after blending this rank's
 * blended data would only cover 0-8s, and so we would not be able
 * to extract the requested output */
static int join_patches(
    int const n_patches, int const *const n_traces, int const *const n_times,
    long int const *const *const shottimes, int const *const *const channels,
    enum AGDTraceType const *const *const trace_types, int const n_patches_out,
    int const *const n_traces_out, int const *const n_times_out,
    long int const *const *const shottimes_out,
    int const *const *const channels_out,
    enum AGDTraceType const *const *const trace_types_out,
    int *const n_patches_joined, int **const n_traces_joined,
    int **const n_times_joined, long int const ***const shottimes_joined,
    int const ***const channels_joined,
    enum AGDTraceType const ***const trace_types_joined) {
  /* Allocate for new joined patches */
  *n_patches_joined = n_patches + n_patches_out;
  *n_traces_joined = (int *)calloc((size_t)*n_patches_joined, sizeof(int *));
  if (*n_traces_joined == NULL) goto err;
  *n_times_joined = (int *)calloc((size_t)*n_patches_joined, sizeof(int *));
  if (*n_times_joined == NULL) goto err;
  *shottimes_joined = (long int const **)calloc((size_t)*n_patches_joined,
                                                sizeof(long int const **));
  if (*shottimes_joined == NULL) goto err;
  *channels_joined =
      (int const **)calloc((size_t)*n_patches_joined, sizeof(int const **));
  if (*channels_joined == NULL) goto err;
  *trace_types_joined = (enum AGDTraceType const **)calloc(
      (size_t)*n_patches_joined, sizeof(enum AGDTraceType const **));
  if (*trace_types_joined == NULL) goto err;

  /* Copy from input and output */
  memcpy(*n_traces_joined, n_traces, (size_t)n_patches * sizeof(int));
  memcpy(*n_traces_joined + n_patches, n_traces_out,
         (size_t)n_patches_out * sizeof(int));
  memcpy(*n_times_joined, n_times, (size_t)n_patches * sizeof(int));
  memcpy(*n_times_joined + n_patches, n_times_out,
         (size_t)n_patches_out * sizeof(int));
  memcpy(*shottimes_joined, shottimes,
         (size_t)n_patches * sizeof(long int const *));
  memcpy(*shottimes_joined + n_patches, shottimes_out,
         (size_t)n_patches_out * sizeof(long int const *));
  memcpy(*channels_joined, channels, (size_t)n_patches * sizeof(int const *));
  memcpy(*channels_joined + n_patches, channels_out,
         (size_t)n_patches_out * sizeof(int const *));
  memcpy(*trace_types_joined, trace_types,
         (size_t)n_patches * sizeof(enum AGDTraceType const *));
  memcpy(*trace_types_joined + n_patches, trace_types_out,
         (size_t)n_patches_out * sizeof(enum AGDTraceType const *));

  return 0;
err:
  fprintf(stderr, "ERROR in join_patches\n");
  return 1;
}
#endif /* AGD_MPI */

int agd_blend(int const n_patches, int const *const n_traces,
              int const *const n_times, long int const *const *const shottimes,
              int const *const *const channels,
              enum AGDTraceType const *const *const trace_types,
              AGD_TYPE *const *const data, enum AGDBlendMode const blend_mode,
              int const taper_length, int const n_patches_out,
              int const *const n_traces_out, int const *const n_times_out,
              long int const *const *const shottimes_out,
              int const *const *const channels_out,
              enum AGDTraceType const *const *const trace_types_out,
#ifdef AGD_MPI
              MPI_Comm comm,
#endif /* AGD_MPI */
              AGD_TYPE *const *const data_out) {
  int patch_idx;
  struct BlendConfig blend_config = ZERO_INIT;
  struct BlendParams *blend_params = NULL;
  struct BlendParams *blend_params_out = NULL;
  AGD_TYPE ***blended = NULL;

  /* Check the inputs */
  if (check_blend_inputs(
          n_patches, n_traces, n_times, shottimes, channels, trace_types, data,
          blend_mode, taper_length, n_patches_out, n_traces_out, n_times_out,
          shottimes_out, channels_out, trace_types_out, data_out))
    goto err;

    /* Set blend config */
#ifdef AGD_MPI
  {
    /* Join the input and output patches before creating blend_config,
     * for the reason described above (before join_patches) */
    int n_patches_joined = 0;
    int *n_traces_joined = NULL;
    int *n_times_joined = NULL;
    long int const **shottimes_joined = NULL;
    int const **channels_joined = NULL;
    enum AGDTraceType const **trace_types_joined = NULL;
    if (join_patches(n_patches, n_traces, n_times, shottimes, channels,
                     trace_types, n_patches_out, n_traces_out, n_times_out,
                     shottimes_out, channels_out, trace_types_out,
                     &n_patches_joined, &n_traces_joined, &n_times_joined,
                     &shottimes_joined, &channels_joined, &trace_types_joined))
      goto err;

    if (set_blend_config(n_patches_joined, n_traces_joined, n_times_joined,
                         shottimes_joined, channels_joined, trace_types_joined,
                         comm, &blend_config))
      goto err;
    free(n_traces_joined);
    free(n_times_joined);
    free(shottimes_joined);
    free(channels_joined);
    free(trace_types_joined);
  }
#else
  if (set_blend_config(n_patches, n_traces, n_times, shottimes, channels,
                       trace_types, &blend_config))
    goto err;
#endif /* AGD_MPI */

  /* Set blend params in */
  blend_params = (struct BlendParams *)malloc((size_t)n_patches *
                                              sizeof(struct BlendParams));
  if (blend_params == NULL) goto err;
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    if (set_blend_params(n_traces[patch_idx], n_times[patch_idx],
                         shottimes[patch_idx], channels[patch_idx],
                         trace_types[patch_idx], &blend_config,
                         blend_params + patch_idx))
      goto err;
  }

  /* Set blend params out */
  blend_params_out = (struct BlendParams *)malloc((size_t)n_patches_out *
                                                  sizeof(struct BlendParams));
  if (blend_params_out == NULL) goto err;
  for (patch_idx = 0; patch_idx < n_patches_out; patch_idx++) {
    if (set_blend_params(n_traces_out[patch_idx], n_times_out[patch_idx],
                         shottimes_out[patch_idx], channels_out[patch_idx],
                         trace_types_out[patch_idx], &blend_config,
                         blend_params_out + patch_idx))
      goto err;
  }

  /* Allocate blended */
  if (allocate_blended(&blend_config, &blended)) goto err;
  zero_blended(&blend_config, (AGD_TYPE *const *const *)blended);

  /* Blend */
  if (blend_mode == AGDBlendSum) {
    for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      blend_sum_forward(data[patch_idx], blend_params + patch_idx,
                        (AGD_TYPE *const *const *)blended);
    }
#ifdef AGD_MPI
    if (blend_sum_forward_mpi(&blend_config, (AGD_TYPE *const *const *)blended))
      goto err;
#endif /* AGD_MPI*/

  } else if (blend_mode == AGDBlendOverwrite) {
    for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      blend_overwrite_forward(data[patch_idx], blend_params + patch_idx,
                              (AGD_TYPE *const *const *)blended);
    }
#ifdef AGD_MPI
    if (blend_overwrite_forward_mpi(&blend_config,
                                    (AGD_TYPE *const *const *)blended))
      goto err;
#endif /* AGD_MPI*/

  } else {
    blend_mean_forward(n_patches, (AGD_TYPE const *const *)data, &blend_config,
                       blend_params, taper_length,
                       (AGD_TYPE *const *const *)blended);
  }

  apply_mute(&blend_config, (AGD_TYPE *const *const *)blended);

  /* Get output */
  for (patch_idx = 0; patch_idx < n_patches_out; patch_idx++) {
    blend_adjoint((AGD_TYPE const *const *const *)blended,
                  blend_params_out + patch_idx, data_out[patch_idx]);
  }

  /* Free allocated memory */
  free_blended(&blend_config, &blended);
  free_blend_params_array(n_patches, &blend_params);
  free_blend_params_array(n_patches_out, &blend_params_out);
  free_blend_config(&blend_config);
  return 0;
err:
  fprintf(stderr, "ERROR in agd_blend\n");
  free_blended(&blend_config, &blended);
  free_blend_params_array(n_patches, &blend_params);
  free_blend_params_array(n_patches_out, &blend_params_out);
  free_blend_config(&blend_config);
  return 1;
}
