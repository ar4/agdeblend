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

/* Set use_taper to 1 for edges that are adjacent to the other patch */
static void set_use_taper_patch(int const patch_volume, int const patch_n_dims,
                                int const *const patch_coords,
                                int const other_volume,
                                int const *const other_coords,
                                int *const patch_use_taper) {
  int dim_idx;

  if (other_volume != patch_volume) return; /* Different volume */

  /* Check if other is a neighbour */
  for (dim_idx = 0; dim_idx < patch_n_dims - 1; dim_idx++) {
    int const patch_coord = patch_coords[dim_idx];
    int const other_coord = other_coords[dim_idx];
    if (other_coord != patch_coord - 1 && other_coord != patch_coord &&
        other_coord != patch_coord + 1) {
      /* Not a neighbour */
      return;
    }
  }

  /* Set use_taper for edges adjacent to other */
  for (dim_idx = 0; dim_idx < patch_n_dims - 1; dim_idx++) {
    int const patch_coord = patch_coords[dim_idx];
    int const other_coord = other_coords[dim_idx];
    if (other_coord == patch_coord - 1) {
      patch_use_taper[dim_idx * 2 + 0] = 1;
    } else if (other_coord == patch_coord + 1) {
      patch_use_taper[dim_idx * 2 + 1] = 1;
    }
  }
}

#ifdef AGD_MPI
/* Set use_taper to 1 for edges adjacent to patches on other ranks */
static int set_use_taper_mpi(int const n_patches, int const *const volumes,
                             int const *const n_dims,
                             int const *const *const coords, MPI_Comm comm,
                             int *const *const use_taper) {
  int *rank_volumes = NULL;
  int **rank_coords = NULL;
  int rank_n_patches = 0;
  int comm_rank;
  int comm_size;
  int rank_idx;
  int n_volumes = get_n_volumes(n_patches, volumes);
  if (MPI_Comm_rank(comm, &comm_rank)) goto err;
  if (MPI_Comm_size(comm, &comm_size)) goto err;
  if (MPI_Allreduce(MPI_IN_PLACE, &n_volumes, 1, MPI_INT, MPI_MAX, comm))
    goto err;

  for (rank_idx = 0; rank_idx < comm_size; rank_idx++) {
    int patch_idx;

    /* Get the number of patches and allocate for these */
    if (rank_idx == comm_rank) {
      rank_n_patches = n_patches;
    }
    if (MPI_Bcast(&rank_n_patches, 1, MPI_INT, rank_idx, comm)) goto err;
    rank_volumes = (int *)malloc((size_t)rank_n_patches * sizeof(int));
    if (rank_volumes == NULL) goto err;
    rank_coords = (int **)calloc((size_t)rank_n_patches, sizeof(int *));
    if (rank_coords == NULL) goto err;

    /* Get the volumes and coords */
    if (rank_idx == comm_rank) {
      memcpy(rank_volumes, volumes, (size_t)n_patches * sizeof(int));
    }
    if (MPI_Bcast(rank_volumes, rank_n_patches, MPI_INT, rank_idx, comm))
      goto err;
    for (patch_idx = 0; patch_idx < rank_n_patches; patch_idx++) {
      rank_coords[patch_idx] = (int *)malloc(
          (size_t)(n_dims[rank_volumes[patch_idx]] - 1) * sizeof(int));
      if (rank_coords[patch_idx] == NULL) goto err;
      if (rank_idx == comm_rank) {
        memcpy(rank_coords[patch_idx], coords[patch_idx],
               (size_t)(n_dims[rank_volumes[patch_idx]] - 1) * sizeof(int));
      }
      if (MPI_Bcast(rank_coords[patch_idx], n_dims[rank_volumes[patch_idx]] - 1,
                    MPI_INT, rank_idx, comm))
        goto err;
    }

    /* Set taper patch */
    if (rank_idx != comm_rank) {
      for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
        int const volume = volumes[patch_idx];
        int const patch_n_dims = n_dims[volume];
        int const *const patch_coords = coords[patch_idx];
        int rank_patch_idx;
        for (rank_patch_idx = 0; rank_patch_idx < rank_n_patches;
             rank_patch_idx++) {
          set_use_taper_patch(
              volume, patch_n_dims, patch_coords, rank_volumes[rank_patch_idx],
              rank_coords[rank_patch_idx], use_taper[patch_idx]);
        }
      }
    }

    /* Free memory */
    free(rank_volumes);
    rank_volumes = NULL;
    if (rank_coords != NULL) {
      for (patch_idx = 0; patch_idx < rank_n_patches; patch_idx++) {
        free(rank_coords[patch_idx]);
        rank_coords[patch_idx] = NULL;
      }
      free(rank_coords);
      rank_coords = NULL;
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_use_taper_mpi\n");
  free(rank_volumes);
  rank_volumes = NULL;
  if (rank_coords != NULL) {
    int patch_idx;
    for (patch_idx = 0; patch_idx < rank_n_patches; patch_idx++) {
      free(rank_coords[patch_idx]);
      rank_coords[patch_idx] = NULL;
    }
    free(rank_coords);
    rank_coords = NULL;
  }
  return 1;
}
#endif /* AGD_MPI */

static int set_use_taper(int const n_patches, int const *const volumes,
                         int const *const n_dims,
                         int const *const *const coords,
#ifdef AGD_MPI
                         MPI_Comm comm,
#endif /* AGD_MPI */
                         int ***const use_taper) {
  int patch_idx;

  *use_taper = (int **)calloc((size_t)n_patches, sizeof(int *));
  if (*use_taper == NULL) goto err;

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const volume = volumes[patch_idx];
    int const patch_n_dims = n_dims[volume];
    int const *const patch_coords = coords[patch_idx];
    int patch_idx2;

    (*use_taper)[patch_idx] =
        (int *)calloc((size_t)(2 * (patch_n_dims - 1)), sizeof(int));
    if ((*use_taper)[patch_idx] == NULL) goto err;

    for (patch_idx2 = 0; patch_idx2 < n_patches; patch_idx2++) {
      if (patch_idx == patch_idx2) continue;
      set_use_taper_patch(volume, patch_n_dims, patch_coords,
                          volumes[patch_idx2], coords[patch_idx2],
                          (*use_taper)[patch_idx]);
    }
  }

#ifdef AGD_MPI
  if (set_use_taper_mpi(n_patches, volumes, n_dims, coords, comm, *use_taper))
    goto err;
#endif /* AGD_MPI */

  return 0;
err:
  fprintf(stderr, "ERROR in set_use_taper\n");
  return 1;
}

static void free_use_taper(int const n_patches, int ***const use_taper) {
  if (*use_taper != NULL) {
    int patch_idx;
    for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      free((*use_taper)[patch_idx]);
      (*use_taper)[patch_idx] = NULL;
    }
    free(*use_taper);
    *use_taper = NULL;
  }
}

/* Extract the length of the time (last) dimension of each patch.
 * If shape_per_patch != 0, the shape argument is assumed to contain an entry
 * for each patch, and is otherwise assumed to only contain an entry for
 * each volume */
static int set_n_times_per_patch(int const n_patches, int const *const volumes,
                                 int const *const n_dims,
                                 int const *const *const shape,
                                 int const shape_per_patch,
                                 int **const n_times) {
  int patch_idx;
  *n_times = (int *)malloc((size_t)n_patches * sizeof(int));
  if (*n_times == NULL) goto err;

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const patch_n_dims = n_dims[volumes[patch_idx]];
    if (shape_per_patch) {
      (*n_times)[patch_idx] = shape[patch_idx][patch_n_dims - 1];
    } else {
      (*n_times)[patch_idx] = shape[volumes[patch_idx]][patch_n_dims - 1];
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_n_times_per_patch\n");
  return 1;
}

/* Given the wavelet_idxs used in each patch, and the length of each wavelet,
 * get the wavelet length for each patch (assumed to be the same for all
 * traces in the same patch) */
static int set_n_times_wavelet_per_patch(
    int const n_patches, int const *const volumes, int const *const n_traces,
    enum AGDTraceType const *const *const trace_types,
    int const *const wavelet_lengths, int const *const *const wavelet_idxs,
    int **const n_times_wavelet_per_patch) {
  int *volume_wavelet_lengths = NULL;
  int patch_idx;

  if (set_volume_wavelet_lengths(n_patches, volumes, n_traces, trace_types,
                                 wavelet_lengths, wavelet_idxs,
                                 &volume_wavelet_lengths))
    goto err;

  *n_times_wavelet_per_patch = (int *)malloc((size_t)n_patches * sizeof(int));
  if (*n_times_wavelet_per_patch == NULL) goto err;

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const volume_idx = volumes[patch_idx];
    (*n_times_wavelet_per_patch)[patch_idx] =
        volume_wavelet_lengths[volume_idx];
  }

  free(volume_wavelet_lengths);
  return 0;
err:
  fprintf(stderr, "ERROR in set_n_times_wavelet_per_patch\n");
  free(volume_wavelet_lengths);
  return 1;
}

static int set_patch_properties(
    int const n_patches, int const *const volumes, int const *const n_dims,
    int const *const *const patch_shapes, int const *const *const window_shapes,
    enum AGDTraceType const *const *const trace_types,
    int const *const wavelet_lengths, int const *const *const wavelet_idxs,
    struct PatchProperties *const patch_properties) {
  /* Set n_traces */
  if (set_n_traces_per_patch(n_patches, volumes, n_dims, patch_shapes,
                             &(patch_properties->n_traces)))
    goto err;

  /* Set trace length */
  if (set_n_times_per_patch(n_patches, volumes, n_dims, patch_shapes, 1,
                            &(patch_properties->n_times_data)))
    goto err;

  /* Set time window length */
  if (set_n_times_per_patch(n_patches, volumes, n_dims, window_shapes, 0,
                            &(patch_properties->n_times_window)))
    goto err;

  /* Set n_times_wavelet */
  if (set_n_times_wavelet_per_patch(
          n_patches, volumes, patch_properties->n_traces, trace_types,
          wavelet_lengths, wavelet_idxs, &(patch_properties->n_times_wavelet)))
    goto err;

  return 0;
err:
  fprintf(stderr, "ERROR in set_patch_properties\n");
  return 1;
}

static void free_patch_properties(
    struct PatchProperties *const patch_properties) {
  free(patch_properties->n_traces);
  patch_properties->n_traces = NULL;
  free(patch_properties->n_times_data);
  patch_properties->n_times_data = NULL;
  free(patch_properties->n_times_window);
  patch_properties->n_times_window = NULL;
  free(patch_properties->n_times_wavelet);
  patch_properties->n_times_wavelet = NULL;
}

/* Set the number of time samples in each patch's traces at different stages
 *
 * All traces in a patch are expected to have the same length as each other
 * at each stage, but traces in different patches may have different lengths.
 *
 * The traces that are output from fk_forward, input to and output from
 * wavelet_forward, and input to blend_forward, have capacity n_times_conv.
 * The length of the data portion of these traces changes, however.
 * The length of the data portion of traces produced by fk_forward is of
 * length n_times_fk, a multiple of the window length in the time
 * dimension.
 * These are convolved with the source wavelets (if provided) to make
 * the output of wavelet_forward, with a data portion of length n_times_conv.
 * n_times_conv (or n_times_fk, if there is no source wavelet) should
 * be at least as long as the input data traces plus
 * the taper length (half a window) at the top and bottom.
 * If wavelet convolution is not used, n_times_fk and n_times_conv are equal.
 */
static int set_n_model_times(int const n_patches, int const *const n_times_data,
                             int const *const n_times_window,
                             int const *const n_times_wavelet,
                             int **const n_times_fk, int **const n_times_conv) {
  int patch_idx;
  *n_times_fk = (int *)malloc((size_t)n_patches * sizeof(int));
  if (*n_times_fk == NULL) goto err;
  *n_times_conv = (int *)malloc((size_t)n_patches * sizeof(int));
  if (*n_times_conv == NULL) goto err;

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const patch_n_times_data = n_times_data[patch_idx];
    int const patch_n_times_window = n_times_window[patch_idx];
    int const patch_n_times_halfwindow = patch_n_times_window / 2;
    int const patch_n_times_wavelet = n_times_wavelet[patch_idx];
    int const n_times_data_before_conv =
        patch_n_times_data - patch_n_times_wavelet + 1;
    int const n_time_windows_without_taper =
        (n_times_data_before_conv + patch_n_times_halfwindow - 1) /
        patch_n_times_halfwindow;
    int const n_time_windows_with_taper = n_time_windows_without_taper + 1;
    (*n_times_fk)[patch_idx] =
        get_combined_time(n_time_windows_with_taper, patch_n_times_window);
    (*n_times_conv)[patch_idx] =
        (*n_times_fk)[patch_idx] + patch_n_times_wavelet - 1;
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_n_model_times\n");
  return 1;
}

static void free_space_window(struct SpaceWindow *const space_window) {
  free(space_window->coords);
  space_window->coords = NULL;
  free(space_window->space_shape);
  space_window->space_shape = NULL;
  free(space_window->taper_length);
  space_window->taper_length = NULL;
  free(space_window->use_taper);
  space_window->use_taper = NULL;
}

static void free_model_window(struct ModelWindow *const window) {
  free_space_window(&(window->space_window));
  free_wavelet_params(&(window->wavelet_params));
  free_blend_params(&(window->blend_params));
}

static void free_model_windows(int const n_windows,
                               struct ModelWindow **const model_windows) {
  int window_idx;

  if (*model_windows != NULL) {
    for (window_idx = 0; window_idx < n_windows; window_idx++) {
      struct ModelWindow *const window = (*model_windows) + window_idx;
      free_model_window(window);
    }
    free(*model_windows);
    *model_windows = NULL;
  }
}

static void free_shottimes(int const n_patches, long int ***const shottimes) {
  if (*shottimes) {
    int patch_idx;
    for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      free((*shottimes)[patch_idx]);
      (*shottimes)[patch_idx] = NULL;
    }
    free(*shottimes);
    *shottimes = NULL;
  }
}

static int allocate_shottimes(int const n_patches, int const *const n_traces,
                              long int ***const shottimes) {
  int patch_idx;

  *shottimes = (long int **)calloc((size_t)n_patches, sizeof(long int *));
  if (*shottimes == NULL) goto err;

  /* Allocate each patch */
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    (*shottimes)[patch_idx] =
        (long int *)malloc((size_t)n_traces[patch_idx] * sizeof(long int));
    if ((*shottimes)[patch_idx] == NULL) goto err;
  }

  return 0;

err:
  fprintf(stderr, "ERROR in allocate_shottimes\n");
  return 1;
}

/* Shift shottimes back to include the time taper */
static int set_shottimes_with_taper(int const n_patches,
                                    int const *const n_traces,
                                    int const *const n_times_window,
                                    long int const *const *const shottimes,
                                    long int ***const shottimes_with_taper) {
  int patch_idx;

  if (allocate_shottimes(n_patches, n_traces, shottimes_with_taper)) goto err;

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const n_times_taper = n_times_window[patch_idx] / 2;
    int trace_idx;
    for (trace_idx = 0; trace_idx < n_traces[patch_idx]; trace_idx++) {
      long int const shottime = shottimes[patch_idx][trace_idx];
      (*shottimes_with_taper)[patch_idx][trace_idx] = shottime - n_times_taper;
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_shottimes_with_taper\n");
  return 1;
}

/* Get the minimum (base) window length: the largest length that we can fit
 * n_windows of (with overlapping) within the data length.
 * Windows are composed of a middle and equal left and right "taper"s.
 * The start positions of adjacent windows will differ by the length of
 * the middle plus one taper length (I imagine it as the left taper).
 * The length of a window is thus this difference in start positions plus
 * a taper length. We need to subtract a taper length from the data length
 * in this calculation to account for the rightmost taper. */
static int get_base_window_length(int const data_length, int const taper_length,
                                  int const n_windows) {
  return (data_length - taper_length) / n_windows + taper_length;
}

/* If the base window length does not cover the whole data length, some of
 * the windows need to be 1 trace larger. This calculates how many. */
static int get_n_larger_windows(int const data_length, int const taper_length,
                                int const base_window_length,
                                int const n_windows) {
  return data_length - taper_length -
         (base_window_length - taper_length) * n_windows;
}

static int get_n_space_windows_in_dim(int const data_length,
                                      int const window_length) {
  if (window_length == data_length) return 1;
  return data_length / (window_length / 2) - 1;
}

static int get_n_space_windows_in_patch(int const n_dims,
                                        int const *const patch_shape,
                                        int const *const window_shape) {
  int dim_idx;
  int n_windows = 1;

  for (dim_idx = 0; dim_idx < n_dims - 1; dim_idx++) {
    int const n_windows_in_dim =
        get_n_space_windows_in_dim(patch_shape[dim_idx], window_shape[dim_idx]);

    n_windows *= n_windows_in_dim;
  }
  return n_windows;
}

static int get_n_space_windows(int const n_patches, int const *const volumes,
                               int const *const n_dims,
                               int const *const *const patch_shapes,
                               int const *const *const window_shapes) {
  int patch_idx;
  int n_windows = 0;
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const patch_n_dims = n_dims[volumes[patch_idx]];
    int const *const patch_patch_shape = patch_shapes[patch_idx];
    int const *const patch_window_shape = window_shapes[volumes[patch_idx]];
    int const n_windows_in_patch = get_n_space_windows_in_patch(
        patch_n_dims, patch_patch_shape, patch_window_shape);
    n_windows += n_windows_in_patch;
  }
  return n_windows;
}

static int allocate_space_window_coordinate_arrays(
    int const n_patches, int const *const volumes, int const *const n_dims,
    int const *const *const patch_shapes, int const *const *const window_shapes,
    struct ModelWindow *const windows) {
  int patch_window_start_idx = 0;
  int patch_idx;
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const patch_n_dims = n_dims[volumes[patch_idx]];
    int const patch_n_space_dims = patch_n_dims - 1;
    int const *const patch_shape = patch_shapes[patch_idx];
    int const *const patch_window_shape = window_shapes[volumes[patch_idx]];
    int const n_windows_in_patch = get_n_space_windows_in_patch(
        patch_n_dims, patch_shape, patch_window_shape);
    int patch_window_idx;
    for (patch_window_idx = 0; patch_window_idx < n_windows_in_patch;
         patch_window_idx++) {
      int const window_idx = patch_window_start_idx + patch_window_idx;
      struct SpaceWindow *const space_window =
          &(windows[window_idx].space_window);
      space_window->coords =
          (int *)malloc((size_t)patch_n_space_dims * sizeof(int));
      if (space_window->coords == NULL) goto err;
      space_window->space_shape =
          (int *)malloc((size_t)patch_n_space_dims * sizeof(int));
      if (space_window->space_shape == NULL) goto err;
      space_window->taper_length =
          (int *)malloc((size_t)patch_n_space_dims * sizeof(int));
      if (space_window->taper_length == NULL) goto err;
      space_window->use_taper =
          (int *)malloc((size_t)(2 * patch_n_space_dims) * sizeof(int));
      if (space_window->use_taper == NULL) goto err;
      space_window->patch_idx = patch_idx;
      space_window->n_dims = patch_n_dims;
    }
    patch_window_start_idx += n_windows_in_patch;
  }

  return 0;
err:
  fprintf(stderr, "ERROR in allocate_space_window_coordinate_arrays\n");
  return 1;
}

static int check_window_shapes_taper(int const n_patches,
                                     int const *const volumes,
                                     int const *const n_dims,
                                     int const *const *const window_shapes,
                                     int const *const *const use_taper) {
  int patch_idx;
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const volume = volumes[patch_idx];
    int const patch_n_dims = n_dims[volume];
    int const *const patch_window_shape = window_shapes[volume];
    int const *const patch_use_taper = use_taper[patch_idx];
    int dim_idx;
    for (dim_idx = 0; dim_idx < patch_n_dims - 1; dim_idx++) {
      int const use_taper_left = patch_use_taper[2 * dim_idx + 0];
      int const use_taper_right = patch_use_taper[2 * dim_idx + 1];
      int const dim_window_shape = patch_window_shape[dim_idx];
      if ((use_taper_left || use_taper_right) &&
          ((dim_window_shape % 2) == 1)) {
        fprintf(stderr,
                "ERROR: window shapes may only be odd if there is no tapering. "
                "Patch %d touches another patch in dimension %d and so will be "
                "tapered, "
                "but has window shape %d.\n",
                patch_idx, dim_idx, dim_window_shape);
        return 1;
      }
    }
  }
  return 0;
}

static void set_space_windows_coordinates(
    int const n_patches, int const *const volumes, int const *const n_dims,
    int const *const *const patch_shapes, int const *const *const window_shapes,
    int const *const *const use_taper, struct ModelWindow *const windows) {
  int patch_window_start_idx = 0;
  int patch_idx;

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const patch_n_dims = n_dims[volumes[patch_idx]];
    int const *const patch_shape = patch_shapes[patch_idx];
    int const *const patch_window_shape = window_shapes[volumes[patch_idx]];
    int const *const patch_use_taper = use_taper[patch_idx];
    int n_windows[AGD_MAX_N_COORDS];
    int n_windows_in_patch = 1;
    int dim_idx;
    for (dim_idx = 0; dim_idx < patch_n_dims - 1; dim_idx++) {
      n_windows[dim_idx] = get_n_space_windows_in_dim(
          patch_shape[dim_idx], patch_window_shape[dim_idx]);
      n_windows_in_patch *= n_windows[dim_idx];
    }
    for (dim_idx = 0; dim_idx < patch_n_dims - 1; dim_idx++) {
      int const taper_length = patch_window_shape[dim_idx] / 2;
      int const n_windows_in_dim = n_windows[dim_idx];
      int const base_window_length = get_base_window_length(
          patch_shape[dim_idx], taper_length, n_windows_in_dim);
      int const n_larger_windows =
          get_n_larger_windows(patch_shape[dim_idx], taper_length,
                               base_window_length, n_windows_in_dim);
      int const n_plane_windows = n_windows_in_patch / n_windows_in_dim;
      int window_start = 0;
      int dim_window_idx;
      for (dim_window_idx = 0; dim_window_idx < n_windows_in_dim;
           dim_window_idx++) {
        int const window_length =
            base_window_length + (dim_window_idx < n_larger_windows);
        int coords[AGD_MAX_N_COORDS] = ZERO_INIT;
        int plane_window_idx;
        coords[dim_idx] = dim_window_idx;

        /* Set for all other windows with the same coordinate in this
         * dimension. */
        for (plane_window_idx = 0; plane_window_idx < n_plane_windows;
             plane_window_idx++) {
          struct SpaceWindow *space_window;
          int window_idx = patch_window_start_idx;
          int dim_stride = 1;
          int coord_dim_idx;

          /* Get the index of the window within the window array */
          for (coord_dim_idx = patch_n_dims - 2; coord_dim_idx >= 0;
               coord_dim_idx--) {
            window_idx += dim_stride * coords[coord_dim_idx];
            dim_stride *= n_windows[coord_dim_idx];
          }

          /* Set the values of the window for this dimension */
          space_window = &(windows[window_idx].space_window);
          space_window->coords[dim_idx] = window_start;
          space_window->space_shape[dim_idx] = window_length;
          space_window->taper_length[dim_idx] = taper_length;
          if (dim_window_idx == 0) {
            space_window->use_taper[2 * dim_idx + 0] =
                patch_use_taper[2 * dim_idx + 0];
          } else {
            space_window->use_taper[2 * dim_idx + 0] = 1;
          }
          if (dim_window_idx == n_windows_in_dim - 1) {
            space_window->use_taper[2 * dim_idx + 1] =
                patch_use_taper[2 * dim_idx + 1];
          } else {
            space_window->use_taper[2 * dim_idx + 1] = 1;
          }

          /* Update the coordinates to the next window in this plane */
          for (coord_dim_idx = patch_n_dims - 2; coord_dim_idx >= 0;
               coord_dim_idx--) {
            if (coord_dim_idx == dim_idx) continue;
            coords[coord_dim_idx]++;
            if (coords[coord_dim_idx] == n_windows[coord_dim_idx]) {
              coords[coord_dim_idx] = 0;
            } else {
              break;
            }
          }
        }
        window_start += window_length - taper_length;
      }
    }
    patch_window_start_idx += n_windows_in_patch;
  }
}

/* Calculate the number of traces in each window, and set it in the struct */
static void set_space_window_n_traces(int const n_windows,
                                      struct ModelWindow *const windows) {
  int window_idx;
  for (window_idx = 0; window_idx < n_windows; window_idx++) {
    struct SpaceWindow *const space_window =
        &(windows[window_idx].space_window);
    int dim_idx;
    space_window->n_traces = 1;
    for (dim_idx = 0; dim_idx < space_window->n_dims - 1; dim_idx++) {
      space_window->n_traces *= space_window->space_shape[dim_idx];
    }
  }
}

/* Divide the space (non-time) shape of patches into overlapping windows of
 * approximately window_shapes shape. */
static int set_space_windows(int const n_patches, int const *const volumes,
                             int const *const n_dims,
                             int const *const *const window_shapes,
                             int const *const *const coords,
                             int const *const *const patch_shapes,
#ifdef AGD_MPI
                             MPI_Comm comm,
#endif /* AGD_MPI */
                             int *const n_windows,
                             struct ModelWindow **const windows) {

  int **use_taper = NULL;

  /* Get the total number of space windows */
  *n_windows = get_n_space_windows(n_patches, volumes, n_dims, patch_shapes,
                                   window_shapes);

  /* Allocate windows */
  *windows = (struct ModelWindow *)calloc((size_t)*n_windows,
                                          sizeof(struct ModelWindow));
  if (*windows == NULL) goto err;

  /* Allocate each space window's coordinate arrays and set properties */
  if (allocate_space_window_coordinate_arrays(
          n_patches, volumes, n_dims, patch_shapes, window_shapes, *windows))
    goto err;

  /* Set use taper */
  if (set_use_taper(n_patches, volumes, n_dims, coords,
#ifdef AGD_MPI
                    comm,
#endif /* AGD_MPI */
                    &use_taper))
    goto err;

  /* Check that window shape is even if using taper */
  if (check_window_shapes_taper(n_patches, volumes, n_dims, window_shapes,
                                (int const *const *)use_taper))
    goto err;

  /* Set each space window's coordinates */
  set_space_windows_coordinates(n_patches, volumes, n_dims, patch_shapes,
                                window_shapes, (int const *const *)use_taper,
                                *windows);

  /* Set the number of traces in each window */
  set_space_window_n_traces(*n_windows, *windows);

  free_use_taper(n_patches, &use_taper);
  return 0;
err:
  fprintf(stderr, "ERROR in set_space_window\n");
  free_use_taper(n_patches, &use_taper);
  return 1;
}

/* Extract a window's portion of the shottimes, channels, etc. arrays */
static int set_window_arrays(struct SpaceWindow const *const window,
                             int const *const *const patch_shapes,
                             long int const *const *const shottimes,
                             int const *const *const channels,
                             enum AGDTraceType const *const *const trace_types,
                             int const *const *const wavelet_idxs,
                             struct WindowArrays *const window_arrays) {
  int const *const coords = window->coords;
  int const *const space_shape = window->space_shape;
  int const patch_idx = window->patch_idx;
  int const n_dims = window->n_dims;
  int const n_traces = window->n_traces;
  int window_trace_idx;

  /* Allocate memory to store the window's portion */
  window_arrays->shottimes =
      (long int *)malloc((size_t)n_traces * sizeof(long int));
  if (window_arrays->shottimes == NULL) goto err;
  window_arrays->channels = (int *)malloc((size_t)n_traces * sizeof(int));
  if (window_arrays->channels == NULL) goto err;
  window_arrays->trace_types =
      (enum AGDTraceType *)malloc((size_t)n_traces * sizeof(enum AGDTraceType));
  if (window_arrays->trace_types == NULL) goto err;
  if (wavelet_idxs != NULL) {
    window_arrays->wavelet_idxs = (int *)malloc((size_t)n_traces * sizeof(int));
    if (window_arrays->wavelet_idxs == NULL) goto err;
  }

  /* Copy into the window's portion */
  for (window_trace_idx = 0; window_trace_idx < n_traces; window_trace_idx++) {
    int patch_trace_idx = 0;
    int dim_stride = 1;
    int window_trace_idx_remainder = window_trace_idx;
    int dim_idx;
    for (dim_idx = n_dims - 2; dim_idx >= 0; dim_idx--) {
      int window_dim_trace_idx =
          window_trace_idx_remainder % space_shape[dim_idx];
      patch_trace_idx += (coords[dim_idx] + window_dim_trace_idx) * dim_stride;
      dim_stride *= patch_shapes[patch_idx][dim_idx];
      window_trace_idx_remainder /= space_shape[dim_idx];
    }
    window_arrays->shottimes[window_trace_idx] =
        shottimes[patch_idx][patch_trace_idx];
    window_arrays->channels[window_trace_idx] =
        channels[patch_idx][patch_trace_idx];
    window_arrays->trace_types[window_trace_idx] =
        trace_types[patch_idx][patch_trace_idx];
    if (wavelet_idxs != NULL) {
      window_arrays->wavelet_idxs[window_trace_idx] =
          wavelet_idxs[patch_idx][patch_trace_idx];
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_window_arrays\n");
  return 1;
}

/* Allocate workspace that is used during forward and backpropagation */
static int allocate_temporary_arrays(struct FKConfigs const *const fk_configs,
                                     void ***const temp_1,
                                     void ***const temp_2) {
  int config_idx;
  size_t n_required_bytes = 0;

  /* Determine the maximum space that will be required */
  for (config_idx = 0; config_idx < fk_configs->n_window_configs;
       config_idx++) {
    struct FKWindowConfig const *const config =
        fk_configs->window_configs + config_idx;
    int n_elements;
    size_t n_bytes;

    /* FK domain */
    n_elements = config->n_time_windows * config->n_elements_padded_window_fk;
    n_bytes = (size_t)n_elements * sizeof(AGD_FFTW_COMPLEX);
    if (n_bytes > n_required_bytes) n_required_bytes = n_bytes;

    /* After FK transformation to TX */
    n_elements = config->n_time_windows * config->n_elements_padded_window;
    n_bytes = (size_t)n_elements * sizeof(AGD_TYPE);
    if (n_bytes > n_required_bytes) n_required_bytes = n_bytes;

    /* After removal of padding and summing of time windows */
    n_elements = config->n_elements_tx_output_capacity;
    n_bytes = (size_t)n_elements * sizeof(AGD_TYPE);
    if (n_bytes > n_required_bytes) n_required_bytes = n_bytes;
  }

  /* Allocate the space */
#ifdef AGD_THREADS
  int const n_threads = omp_get_max_threads();
#else
  int const n_threads = 1;
#endif /* AGD_THREADS */
  int thread_idx;
  *temp_1 = (void **)calloc((size_t)n_threads, sizeof(void *));
  if (*temp_1 == NULL) goto err;
  *temp_2 = (void **)calloc((size_t)n_threads, sizeof(void *));
  if (*temp_2 == NULL) goto err;
  if (n_required_bytes > 0) {
    for (thread_idx = 0; thread_idx < n_threads; thread_idx++) {
      (*temp_1)[thread_idx] = (void *)AGD_FFTW_MALLOC(n_required_bytes);
      if ((*temp_1)[thread_idx] == NULL) goto err;
      (*temp_2)[thread_idx] = (void *)AGD_FFTW_MALLOC(n_required_bytes);
      if ((*temp_2)[thread_idx] == NULL) goto err;
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in allocate_temporary_arrays\n");
  return 1;
}

static void free_temporary_arrays(void ***const temp_1, void ***const temp_2) {
#ifdef AGD_THREADS
  int const n_threads = omp_get_max_threads();
#else
  int const n_threads = 1;
#endif /* AGD_THREADS */
  int thread_idx;
  for (thread_idx = 0; thread_idx < n_threads; thread_idx++) {
    AGD_FFTW_FREE((*temp_1)[thread_idx]);
    (*temp_1)[thread_idx] = NULL;
    AGD_FFTW_FREE((*temp_2)[thread_idx]);
    (*temp_2)[thread_idx] = NULL;
  }
  free(*temp_1);
  *temp_1 = NULL;
  free(*temp_2);
  *temp_2 = NULL;
}

static void free_window_arrays(struct WindowArrays *const window_arrays) {
  free(window_arrays->shottimes);
  window_arrays->shottimes = NULL;
  free(window_arrays->channels);
  window_arrays->channels = NULL;
  free(window_arrays->trace_types);
  window_arrays->trace_types = NULL;
  free(window_arrays->wavelet_idxs);
  window_arrays->wavelet_idxs = NULL;
}

static void free_model(struct Model *const model) {
  free(model->n_times_fk);
  model->n_times_fk = NULL;
  free(model->n_times_conv);
  model->n_times_conv = NULL;
  free_fk_config(&(model->fk_configs));
  free_blend_config(&(model->blend_config));
  free_wavelet_config(&(model->wavelet_config));
  free_model_windows(model->n_windows, &(model->windows));
  model->n_windows = 0;
  free_temporary_arrays(&(model->temp_1), &(model->temp_2));
}

static int contains_live_traces(int const n_traces,
                                enum AGDTraceType const *const trace_types) {
  int trace_idx;
  for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
    if (trace_types[trace_idx] == AGDLive) return 1;
  }
  return 0;
}

static void delete_window(int const delete_window_idx, int *const n_windows,
                          struct ModelWindow *const windows) {
  int window_idx;
  free_model_window(windows + delete_window_idx);
  for (window_idx = delete_window_idx; window_idx < *n_windows - 1;
       window_idx++) {
    windows[window_idx] = windows[window_idx + 1];
  }
  (*n_windows)--;
}

static int set_model(
    int const n_patches, int const *const volumes, int const *const n_dims,
    int const *const *const window_shapes, int const *const *const coords,
    int const *const *const patch_shapes,
    long int const *const *const shottimes, int const *const *const channels,
    enum AGDTraceType const *const *const trace_types,
    int const *const wavelet_lengths, int const *const *const wavelet_idxs,
    AGD_TYPE const *const *const wavelets,
#ifdef AGD_MPI
    MPI_Comm comm,
#endif /* AGD_MPI */
    struct Model *const model) {
  int window_idx;
  struct PatchProperties patch_properties = ZERO_INIT;
  long int **shottimes_with_taper = NULL;

  if (set_patch_properties(n_patches, volumes, n_dims, patch_shapes,
                           window_shapes, trace_types, wavelet_lengths,
                           wavelet_idxs, &patch_properties))
    goto err;

  /* Set n_times_fk, n_times_conv */
  if (set_n_model_times(n_patches, patch_properties.n_times_data,
                        patch_properties.n_times_window,
                        patch_properties.n_times_wavelet, &(model->n_times_fk),
                        &(model->n_times_conv)))
    goto err;

  /* Add taper to shottimes */
  if (set_shottimes_with_taper(n_patches, patch_properties.n_traces,
                               patch_properties.n_times_window, shottimes,
                               &shottimes_with_taper))
    goto err;

  /* Set blend_config */
  if (set_blend_config(
          n_patches, patch_properties.n_traces, model->n_times_conv,
          (long int const *const *)shottimes_with_taper, channels, trace_types,
#ifdef AGD_MPI
          comm,
#endif /* AGD_MPI */
          &(model->blend_config)))
    goto err;

  /* Set space windows */
  set_space_windows(n_patches, volumes, n_dims, window_shapes, coords,
                    patch_shapes,
#ifdef AGD_MPI
                    comm,
#endif /* AGD_MPI */
                    &(model->n_windows), &(model->windows));

  /* Loop over windows */
  for (window_idx = 0; window_idx < model->n_windows; window_idx++) {
    struct WindowArrays window_arrays = ZERO_INIT;
    struct ModelWindow *const window = model->windows + window_idx;
    int const patch_idx = window->space_window.patch_idx;

    /* Set arrays for window */
    if (set_window_arrays(&(window->space_window), patch_shapes,
                          (long int const *const *)shottimes_with_taper,
                          channels, trace_types, wavelet_idxs,
                          &window_arrays)) {
      free_window_arrays(&window_arrays);
      goto err;
    }

    /* Delete the window if it does not contain live traces */
    if (!contains_live_traces(window->space_window.n_traces,
                              window_arrays.trace_types)) {
      free_window_arrays(&window_arrays);
      delete_window(window_idx, &(model->n_windows), model->windows);
      window_idx--;
      continue;
    }

    /* Setup each of the operators used in the model */
    if (set_fk_configs(
            model->n_times_fk[patch_idx], model->n_times_conv[patch_idx],
            patch_properties.n_times_window[patch_idx],
            window->space_window.n_dims, window->space_window.space_shape,
            window->space_window.use_taper, &(model->fk_configs),
            &(window->fk_config_idx))) {
      free_window_arrays(&window_arrays);
      goto err;
    }

    if (set_wavelet_params(
            window->space_window.n_traces, model->n_times_conv[patch_idx],
            window_arrays.trace_types, window_arrays.wavelet_idxs,
            &(model->wavelet_config), &(window->wavelet_params))) {
      free_window_arrays(&window_arrays);
      goto err;
    }

    if (set_blend_params(window->space_window.n_traces,
                         model->n_times_conv[patch_idx],
                         window_arrays.shottimes, window_arrays.channels,
                         window_arrays.trace_types, &(model->blend_config),
                         &(window->blend_params))) {
      free_window_arrays(&window_arrays);
      goto err;
    }

    free_window_arrays(&window_arrays);
  }

  /* Determine the required size of the temporary arrays and allocate them */
  if (allocate_temporary_arrays(&(model->fk_configs), &(model->temp_1),
                                &(model->temp_2)))
    goto err;

  /* Create FFTW plans */
  if (set_fk_fftw_plans(model->temp_1[0], model->temp_2[0],
                        &(model->fk_configs)))
    goto err;
  if (set_wavelet_fftw_plans(model->temp_1[0], model->temp_2[0],
                             &(model->wavelet_config)))
    goto err;

  /* Transform wavelets (if used) */
  if (transform_wavelets(&(model->wavelet_config), wavelet_lengths, wavelets,
                         model->temp_1[0]))
    goto err;

  free_shottimes(n_patches, &shottimes_with_taper);
  free_patch_properties(&patch_properties);
  return 0;
err:
  fprintf(stderr, "ERROR in set_model\n");
  free_shottimes(n_patches, &shottimes_with_taper);
  free_patch_properties(&patch_properties);
  return 1;
}

static int model_forward(AGD_FFTW_COMPLEX const *const *const model_fk,
                         struct Model *const model,
                         AGD_TYPE *const *const *const blended) {
  zero_blended(&(model->blend_config), blended);

#ifdef AGD_THREADS
#pragma omp parallel
  {
    int const thread_idx = omp_get_thread_num();
    void *const temp_1 = model->temp_1[thread_idx];
    void *const temp_2 = model->temp_2[thread_idx];
    int window_idx;
#pragma omp for
#else
  {
    void *const temp_1 = model->temp_1[0];
    void *const temp_2 = model->temp_2[0];
    int window_idx;
#endif /* AGD_THREADS */
    for (window_idx = 0; window_idx < model->n_windows; window_idx++) {
      struct ModelWindow const *const window = model->windows + window_idx;
      struct FKWindowConfig const *const fk_config =
          model->fk_configs.window_configs + window->fk_config_idx;
      AGD_FFTW_COMPLEX const *const window_model_fk = model_fk[window_idx];

      fk_forward(window_model_fk, fk_config, temp_2, temp_1);

      wavelet_forward(&(model->wavelet_config), &(window->wavelet_params),
                      temp_1, temp_2);

#ifdef AGD_THREADS
#pragma omp critical
#endif /* AGD_THREADS */
      blend_sum_forward((AGD_TYPE *)temp_1, &(window->blend_params), blended);
    }
  }
#ifdef AGD_MPI
  if (blend_sum_forward_mpi(&(model->blend_config), blended)) goto err;
#endif /* AGD_MPI*/
  apply_mute(&(model->blend_config), blended);

  return 0;
#ifdef AGD_MPI
err:
  fprintf(stderr, "ERROR in model_forward\n");
  return 1;
#endif /* AGD_MPI*/
}

static void model_adjoint_update(AGD_TYPE const *const *const *const blended,
                                 struct Model *const model,
                                 AGD_TYPE const stepsize,
                                 AGD_FFTW_COMPLEX *const *const model_fk) {
#ifdef AGD_THREADS
#pragma omp parallel
  {
    int const thread_idx = omp_get_thread_num();
    void *const temp_1 = model->temp_1[thread_idx];
    void *const temp_2 = model->temp_2[thread_idx];
    int window_idx;
#pragma omp for
#else
  {
    void *const temp_1 = model->temp_1[0];
    void *const temp_2 = model->temp_2[0];
    int window_idx;
#endif /* AGD_THREADS */
    for (window_idx = 0; window_idx < model->n_windows; window_idx++) {
      struct ModelWindow const *const window = model->windows + window_idx;
      struct FKWindowConfig const *const fk_config =
          model->fk_configs.window_configs + window->fk_config_idx;
      AGD_FFTW_COMPLEX *const window_model_fk = model_fk[window_idx];
      int idx;

      blend_adjoint(blended, &(window->blend_params), (AGD_TYPE *)temp_1);

      wavelet_adjoint(&(model->wavelet_config), &(window->wavelet_params),
                      temp_1, temp_2);

      fk_adjoint(fk_config, temp_2, temp_1);

      /* Update fk_data */
      for (idx = 0; idx < 2 * fk_config->n_time_windows *
                              fk_config->n_elements_padded_window_fk;
           idx++) {
        ((AGD_TYPE *)window_model_fk)[idx] -=
            stepsize * ((AGD_TYPE *)(temp_1))[idx];
      }
    }
  }
}

/* Allocate and zero space to store the model values (in the FK domain) */
static int allocate_model_fk(struct Model const *const model,
                             AGD_FFTW_COMPLEX ***const model_fk) {
  int window_idx;
  *model_fk = (AGD_FFTW_COMPLEX **)calloc((size_t)model->n_windows,
                                          sizeof(AGD_FFTW_COMPLEX *));
  if (*model_fk == NULL) goto err;

  for (window_idx = 0; window_idx < model->n_windows; window_idx++) {
    struct ModelWindow const *const window = model->windows + window_idx;
    struct FKWindowConfig const *const fk_config =
        model->fk_configs.window_configs + window->fk_config_idx;
    size_t const n_elements = (size_t)(fk_config->n_time_windows *
                                       fk_config->n_elements_padded_window_fk);

    (*model_fk)[window_idx] =
        (AGD_FFTW_COMPLEX *)AGD_FFTW_ALLOC_COMPLEX(n_elements);
    if ((*model_fk)[window_idx] == NULL) goto err;

    memset((*model_fk)[window_idx], 0, n_elements * sizeof(AGD_FFTW_COMPLEX));
  }

  return 0;
err:
  fprintf(stderr, "ERROR in allocate_model_fk\n");
  return 1;
}

static void free_model_fk(int const n_windows,
                          AGD_FFTW_COMPLEX ***const model_fk) {
  int window_idx;

  if (*model_fk != NULL) {
    for (window_idx = 0; window_idx < n_windows; window_idx++) {
      AGD_FFTW_FREE((*model_fk)[window_idx]);
      (*model_fk)[window_idx] = NULL;
    }
    free(*model_fk);
    *model_fk = NULL;
  }
}

/* Find the maximum amplitude in the model, as the starting threshold */
static AGD_TYPE get_initial_lamb(
    struct Model const *const model,
    AGD_FFTW_COMPLEX const *const *const model_fk) {
  AGD_TYPE lamb_sq = AGD_ZERO;
  int window_idx;
  for (window_idx = 0; window_idx < model->n_windows; window_idx++) {
    struct ModelWindow const *const window = model->windows + window_idx;
    struct FKWindowConfig const *const fk_config =
        model->fk_configs.window_configs + window->fk_config_idx;
    AGD_FFTW_COMPLEX const *const window_model_fk = model_fk[window_idx];
    int idx;
    for (idx = 0; idx < fk_config->n_time_windows *
                            fk_config->n_elements_padded_window_fk;
         idx++) {
      AGD_TYPE const abs_sq =
          window_model_fk[idx][0] * window_model_fk[idx][0] +
          window_model_fk[idx][1] * window_model_fk[idx][1];
      if (abs_sq > lamb_sq) lamb_sq = abs_sq;
    }
  }
#ifdef AGD_MPI
  MPI_Allreduce(MPI_IN_PLACE, &lamb_sq, 1, AGD_MPI_TYPE, MPI_MAX,
                model->blend_config.comm);
#endif /* AGD_MPI*/
  return AGD_SQRT(lamb_sq);
}

/* Return half the reciprocal of the maximum blending factor as the stepsize */
static AGD_TYPE get_stepsize(struct Model *const model,
                             AGD_TYPE *const *const *const model_blended) {
  struct BlendConfig const *const blend_config = &(model->blend_config);
  AGD_TYPE blend_factor = AGD_ONE;
  int channel_idx;
  int window_idx;

  /* Blend-sum an array of all ones for each window */
  zero_blended(blend_config, model_blended);
  for (window_idx = 0; window_idx < model->n_windows; window_idx++) {
    struct ModelWindow const *const window = model->windows + window_idx;
    struct SpaceWindow const *const space_window = &(window->space_window);
    struct BlendParams const *const blend_params = &(window->blend_params);
    struct FKWindowConfig const *const fk_config =
        model->fk_configs.window_configs + window->fk_config_idx;
    int const patch_idx = space_window->patch_idx;
    int const n_times = model->n_times_conv[patch_idx];
    int const n_traces = space_window->n_traces;
    int const n_dims = space_window->n_dims;
    AGD_TYPE *const temp = (AGD_TYPE *)model->temp_1[0];
    int coords[AGD_MAX_N_COORDS] = ZERO_INIT;
    int trace_idx;
    for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
      AGD_TYPE scaler = AGD_ONE;
      int dim_idx;
      int time_idx;
      AGD_TYPE *const trace = temp + trace_idx * n_times;

      for (dim_idx = n_dims - 2; dim_idx >= 0; dim_idx--) {
        scaler *= fk_config->taper_windows[dim_idx][coords[dim_idx]];
      }

      for (time_idx = 0; time_idx < n_times; time_idx++) {
        trace[time_idx] = scaler;
      }

      for (dim_idx = n_dims - 2; dim_idx >= 0; dim_idx--) {
        coords[dim_idx]++;
        if (coords[dim_idx] == fk_config->window_shape[dim_idx]) {
          coords[dim_idx] = 0;
        } else {
          break;
        }
      }
    }
    blend_sum_forward(temp, blend_params, model_blended);
  }
#ifdef AGD_MPI
  if (blend_sum_forward_mpi(&(model->blend_config), model_blended)) {
    fprintf(stderr, "ERROR reported but ignored in get_stepsize\n");
  }
#endif /* AGD_MPI*/

  /* Find the maximum value in the resulting blended data, this is the blend
   * factor */
  for (channel_idx = 0; channel_idx < blend_config->n_channels; channel_idx++) {
    struct ChannelIntervals const *const channel_intervals =
        blend_config->channels_intervals + channel_idx;
    AGD_TYPE *const *const blended_channel = model_blended[channel_idx];
    int interval_idx;
    for (interval_idx = 0; interval_idx < channel_intervals->n_intervals;
         interval_idx++) {
      struct Interval const *const interval =
          channel_intervals->intervals + interval_idx;
      AGD_TYPE *const blended_interval = blended_channel[interval_idx];
      size_t const n_times = (size_t)(interval->stop - interval->start);
      size_t time_idx;
      for (time_idx = 0; time_idx < n_times; time_idx++) {
        if (blended_interval[time_idx] > blend_factor) {
          blend_factor = blended_interval[time_idx];
        }
      }
    }
  }

#ifdef AGD_MPI
  MPI_Allreduce(MPI_IN_PLACE, &blend_factor, 1, AGD_MPI_TYPE, MPI_MAX,
                blend_config->comm);
#endif /* AGD_MPI*/

  /* Set the stepsize as half the reciprocal of the blend factor */
  return AGD_HALF / blend_factor;
}

static void set_residual(AGD_TYPE const *const *const *const blended,
                         struct BlendConfig const *const blend_config,
                         AGD_TYPE *const *const *const model_blended) {
  int channel_idx;
#ifdef AGD_THREADS
#pragma omp parallel for
#endif /* AGD_THREADS */
  for (channel_idx = 0; channel_idx < blend_config->n_channels; channel_idx++) {
    struct ChannelIntervals const *const channel_intervals =
        blend_config->channels_intervals + channel_idx;
    AGD_TYPE const *const *const blended_channel = blended[channel_idx];
    AGD_TYPE *const *const model_blended_channel = model_blended[channel_idx];
    int interval_idx;
    for (interval_idx = 0; interval_idx < channel_intervals->n_intervals;
         interval_idx++) {
      struct Interval const *const interval =
          channel_intervals->intervals + interval_idx;
      AGD_TYPE const *const blended_interval = blended_channel[interval_idx];
      AGD_TYPE *const model_blended_interval =
          model_blended_channel[interval_idx];
      size_t const n_times = (size_t)(interval->stop - interval->start);
      size_t time_idx;
      for (time_idx = 0; time_idx < n_times; time_idx++) {
        model_blended_interval[time_idx] -= blended_interval[time_idx];
      }
    }
  }
}

static AGD_TYPE get_residual_norm(
    struct BlendConfig const *const blend_config,
    AGD_TYPE const *const *const *const model_blended) {
  int channel_idx;
  AGD_TYPE residual = AGD_ZERO;
#ifdef AGD_THREADS
#pragma omp parallel for reduction(+ : residual)
#endif /* AGD_THREADS */
  for (channel_idx = 0; channel_idx < blend_config->n_channels; channel_idx++) {
    struct ChannelIntervals const *const channel_intervals =
        blend_config->channels_intervals + channel_idx;
    AGD_TYPE const *const *const model_blended_channel =
        model_blended[channel_idx];
    int interval_idx;
    for (interval_idx = 0; interval_idx < channel_intervals->n_intervals;
         interval_idx++) {
      struct Interval const *const interval =
          channel_intervals->intervals + interval_idx;
      AGD_TYPE const *const model_blended_interval =
          model_blended_channel[interval_idx];
      size_t const n_times = (size_t)(interval->stop - interval->start);
      size_t time_idx;
      for (time_idx = 0; time_idx < n_times; time_idx++) {
        residual +=
            model_blended_interval[time_idx] * model_blended_interval[time_idx];
      }
    }
  }
#ifdef AGD_MPI
  {
    MPI_Comm comm = blend_config->comm;
    int comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    if (comm_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, &residual, 1, AGD_MPI_TYPE, MPI_SUM, 0, comm);
    } else {
      MPI_Reduce(&residual, NULL, 1, AGD_MPI_TYPE, MPI_SUM, 0, comm);
    }
  }
#endif /* AGD_MPI */
  return AGD_SQRT(residual);
}

/*
 * if |d| < lamb:
 *   d = 0
 * else:
 *   d = d/|d|(|d| - lamb)
 */
static void shrink_threshold(struct Model const *const model,
                             AGD_TYPE const lamb,
                             AGD_FFTW_COMPLEX *const *const model_fk) {
  AGD_TYPE const lamb_sq = lamb * lamb;
  int window_idx;
#ifdef AGD_THREADS
#pragma omp parallel for
#endif /* AGD_THREADS */
  for (window_idx = 0; window_idx < model->n_windows; window_idx++) {
    struct ModelWindow const *const window = model->windows + window_idx;
    struct FKWindowConfig const *const fk_config =
        model->fk_configs.window_configs + window->fk_config_idx;
    AGD_FFTW_COMPLEX *const window_model_fk = model_fk[window_idx];
    int idx;
    for (idx = 0; idx < fk_config->n_time_windows *
                            fk_config->n_elements_padded_window_fk;
         idx++) {
      AGD_TYPE const abs_sq =
          window_model_fk[idx][0] * window_model_fk[idx][0] +
          window_model_fk[idx][1] * window_model_fk[idx][1];
      if (abs_sq <= lamb_sq) {
        window_model_fk[idx][0] = AGD_ZERO;
        window_model_fk[idx][1] = AGD_ZERO;
      } else {
        AGD_TYPE const abs = AGD_SQRT(abs_sq);
        AGD_TYPE const shrink = (abs - lamb) / abs;
        window_model_fk[idx][0] *= shrink;
        window_model_fk[idx][1] *= shrink;
      }
    }
  }
}

static void set_output(AGD_FFTW_COMPLEX const *const *const model_fk,
                       int const n_patches, int const *const volumes,
                       int const *const n_dims,
                       int const *const *const patch_shapes,
                       struct Model *const model, AGD_TYPE *const *const data) {
  /* Zero output */
  {
    int patch_idx;
    for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      size_t n_elements = 1;
      int dim_idx;
      for (dim_idx = 0; dim_idx < n_dims[volumes[patch_idx]]; dim_idx++) {
        n_elements *= (size_t)patch_shapes[patch_idx][dim_idx];
      }
      memset(data[patch_idx], 0, n_elements * sizeof(AGD_TYPE));
    }
  }

#ifdef AGD_THREADS
#pragma omp parallel
  {
    int const thread_idx = omp_get_thread_num();
    void *const temp_1 = model->temp_1[thread_idx];
    void *const temp_2 = model->temp_2[thread_idx];
    int window_idx;
#pragma omp for
#else
  {
    void *const temp_1 = model->temp_1[0];
    void *const temp_2 = model->temp_2[0];
    int window_idx;
#endif /* AGD_THREADS */
       /* Add windows to output */
    for (window_idx = 0; window_idx < model->n_windows; window_idx++) {
      struct ModelWindow const *const window = model->windows + window_idx;
      struct FKWindowConfig const *const fk_config =
          model->fk_configs.window_configs + window->fk_config_idx;
      struct SpaceWindow const *const space_window = &(window->space_window);
      int const patch_idx = space_window->patch_idx;
      int const *const space_shape = space_window->space_shape;
      int const *const coords = space_window->coords;
      int const *const patch_patch_shape = patch_shapes[patch_idx];
      AGD_FFTW_COMPLEX const *const window_model_fk = model_fk[window_idx];
      AGD_TYPE *const patch_data = data[patch_idx];
      int const patch_n_dims = space_window->n_dims;
      int const n_traces = space_window->n_traces;
      int const patch_n_times = patch_patch_shape[patch_n_dims - 1];
      int const window_n_times = fk_config->n_times_output_capacity;
      int const taper_len = fk_config->window_shape[patch_n_dims - 1] / 2;
      int window_trace_idx;

      fk_forward(window_model_fk, fk_config, temp_2, temp_1);

      wavelet_forward(&(model->wavelet_config), &(window->wavelet_params),
                      temp_1, temp_2);

#ifdef AGD_THREADS
#pragma omp critical
#endif /* AGD_THREADS */
      for (window_trace_idx = 0; window_trace_idx < n_traces;
           window_trace_idx++) {
        int patch_trace_idx = 0;
        int dim_stride = 1;
        int window_trace_idx_remainder = window_trace_idx;
        AGD_TYPE *data_trace;
        AGD_TYPE *window_trace =
            (AGD_TYPE *)temp_1 +
            (size_t)window_trace_idx * (size_t)window_n_times +
            (size_t)taper_len;
        int dim_idx;
        int time_idx;

        /* Determine the position of the trace in the output */
        for (dim_idx = patch_n_dims - 2; dim_idx >= 0; dim_idx--) {
          int window_dim_trace_idx =
              window_trace_idx_remainder % space_shape[dim_idx];
          patch_trace_idx +=
              (coords[dim_idx] + window_dim_trace_idx) * dim_stride;
          dim_stride *= patch_patch_shape[dim_idx];
          window_trace_idx_remainder /= space_shape[dim_idx];
        }
        data_trace =
            patch_data + (size_t)patch_trace_idx * (size_t)patch_n_times;

        for (time_idx = 0; time_idx < patch_n_times; time_idx++) {
          data_trace[time_idx] += window_trace[time_idx];
        }
      }
    }
  }
}
