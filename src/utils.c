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

/* Get the number of traces in one patch */
static int get_n_traces(int const n_dims, int const *const data_shape) {
  int dim_idx;
  int n_traces = 0;

  for (dim_idx = 0; dim_idx < n_dims - 1; dim_idx++) {
    int const dim_shape = data_shape[dim_idx];
    if (n_traces == 0) {
      n_traces = dim_shape;
    } else {
      n_traces *= dim_shape;
    }
  }

  return n_traces;
}

/* Get the number of traces in each patch */
static int set_n_traces_per_patch(int const n_patches, int const *const volumes,
                                  int const *const n_dims,
                                  int const *const *const data_shapes,
                                  int **n_traces_per_patch) {
  int patch_idx;

  *n_traces_per_patch = (int *)malloc((size_t)n_patches * sizeof(int));
  if (*n_traces_per_patch == NULL) goto err;

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    (*n_traces_per_patch)[patch_idx] =
        get_n_traces(n_dims[volumes[patch_idx]], data_shapes[patch_idx]);
  }

  return 0;

err:
  fprintf(stderr, "ERROR in set_n_traces_per_patch\n");
  return 1;
}

/* Get the number of volumes (assumed to be the largest value in volumes + 1) */
static int get_n_volumes(int const n_patches, int const *const volumes) {
  int n_volumes = 0;
  int patch_idx;
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    if (volumes[patch_idx] > n_volumes) n_volumes = volumes[patch_idx];
  }

  return n_volumes + 1;
}

/* Find the wavelet length for each volume.
 * It is assumed to be the same for every trace in the volume, so just use
 * the value from the first non-missing trace */
static int set_volume_wavelet_lengths(
    int const n_patches, int const *const volumes, int const *const n_traces,
    enum AGDTraceType const *const *const trace_types,
    int const *const wavelet_lengths, int const *const *const wavelet_idxs,
    int **const volume_wavelet_lengths) {
  int const n_volumes = get_n_volumes(n_patches, volumes);

  *volume_wavelet_lengths = (int *)calloc((size_t)n_volumes, sizeof(int));
  if (*volume_wavelet_lengths == NULL) goto err;

  if (wavelet_idxs == NULL) {
    /* Not using wavelets - equivalent to convolution with a delta */
    int volume_idx;
    for (volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      (*volume_wavelet_lengths)[volume_idx] = 1;
    }
  } else {
    int volume_idx;
    int patch_idx;
    for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      volume_idx = volumes[patch_idx];
      if ((*volume_wavelet_lengths)[volume_idx] == 0) {
        /* The length for this volume is not yet set.
         * Find the first live trace to get a volume wavelet length */
        int trace_idx;
        for (trace_idx = 0; trace_idx < n_traces[patch_idx]; trace_idx++) {
          if (trace_types[patch_idx][trace_idx] == AGDLive) {
            (*volume_wavelet_lengths)[volume_idx] =
                wavelet_lengths[wavelet_idxs[patch_idx][trace_idx]];
            if ((*volume_wavelet_lengths)[volume_idx] <= 0) {
              fprintf(stderr, "ERROR: Wavelet lengths must be positive\n");
              goto err;
            }
            break;
          }
        }
      }
    }
    /* Set the volume wavelet length to 1 for any volumes that are still unset
     * (they do not have any live traces) */
    for (volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      if ((*volume_wavelet_lengths)[volume_idx] == 0) {
        (*volume_wavelet_lengths)[volume_idx] = 1;
      }
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_volume_wavelet_lengths\n");
  return 1;
}
