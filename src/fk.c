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

static int get_combined_time(int const n_time_windows,
                             int const window_length) {
  return (n_time_windows + 1) * (window_length / 2);
}

static int get_n_time_windows(int const n_times, int const window_length) {
  return n_times / (window_length / 2) - 1;
}

static int set_taper_window(int const window_length, int const *const use_taper,
                            AGD_TYPE **const taper_window) {
  int const taper_length = window_length / 2;
  int idx;

  *taper_window = (AGD_TYPE *)malloc((size_t)window_length * sizeof(AGD_TYPE));
  if (*taper_window == NULL) goto err;

  /* Taper at start */
  if (use_taper[0]) {
    for (idx = 0; idx < taper_length; idx++) {
      (*taper_window)[idx] =
          AGD_HALF *
          (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI * (AGD_TYPE)idx /
                             (AGD_TYPE)(2 * taper_length)));
    }
  } else {
    for (idx = 0; idx < taper_length; idx++) {
      (*taper_window)[idx] = AGD_ONE;
    }
  }

  /* Middle */
  for (idx = taper_length; idx < window_length - taper_length; idx++) {
    (*taper_window)[idx] = AGD_ONE;
  }

  /* Taper at end */
  if (use_taper[1]) {
    for (idx = 0; idx < taper_length; idx++) {
      (*taper_window)[window_length - taper_length + idx] =
          AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                        (AGD_TYPE)(taper_length + idx) /
                                        (AGD_TYPE)(2 * taper_length)));
    }
  } else {
    for (idx = window_length - taper_length; idx < window_length; idx++) {
      (*taper_window)[idx] = AGD_ONE;
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_taper_window\n");
  return 1;
}

static int set_fk_window_config(int const n_times,
                                int const n_times_output_capacity,
                                int const n_times_window, int const n_dims,
                                int const *const window_shape,
                                int const *const use_taper,
                                struct FKWindowConfig *const config) {
  int dim_idx;

  /* n_time_windows */
  config->n_time_windows = get_n_time_windows(n_times, n_times_window);

  /* n_times */
  config->n_times = n_times;

  /* n_times_output_capacity */
  config->n_times_output_capacity = n_times_output_capacity;

  /* n_dims */
  config->n_dims = n_dims;

  /* window_shape */
  config->window_shape = (int *)malloc((size_t)n_dims * sizeof(int));
  if (config->window_shape == NULL) goto err;
  memcpy(config->window_shape, window_shape,
         (size_t)(n_dims - 1) * sizeof(int));
  config->window_shape[n_dims - 1] = n_times_window;

  /* n_traces, n_elements_tx, n_elements_tx_output_capacity */
  config->n_traces = 1;
  for (dim_idx = 0; dim_idx < n_dims - 1; dim_idx++) {
    config->n_traces *= config->window_shape[dim_idx];
  }
  config->n_elements_tx = config->n_traces * config->n_times;
  config->n_elements_tx_output_capacity =
      config->n_traces * config->n_times_output_capacity;

  /* padded_window_shape */
  config->padded_window_shape = (int *)malloc((size_t)n_dims * sizeof(int));
  if (config->padded_window_shape == NULL) goto err;
  for (dim_idx = 0; dim_idx < n_dims; dim_idx++) {
    config->padded_window_shape[dim_idx] = (int)AGD_CEIL(
        (AGD_TYPE)config->window_shape[dim_idx] * (1 + AGD_PAD_FRACTION));
  }

  /* n_traces_padded, n_elements_padded_window */
  config->n_traces_padded = 1;
  for (dim_idx = 0; dim_idx < n_dims - 1; dim_idx++) {
    config->n_traces_padded *= config->padded_window_shape[dim_idx];
  }
  config->n_elements_padded_window =
      config->n_traces_padded * config->padded_window_shape[n_dims - 1];

  /* padded_window_shape_fk */
  config->padded_window_shape_fk = (int *)malloc((size_t)n_dims * sizeof(int));
  if (config->padded_window_shape_fk == NULL) goto err;
  memcpy(config->padded_window_shape_fk, config->padded_window_shape,
         (size_t)(n_dims - 1) * sizeof(int));
  config->padded_window_shape_fk[n_dims - 1] =
      config->padded_window_shape[n_dims - 1] / 2 + 1;

  /* n_elements_padded_window_fk */
  config->n_elements_padded_window_fk = 1;
  for (dim_idx = 0; dim_idx < n_dims; dim_idx++) {
    config->n_elements_padded_window_fk *=
        config->padded_window_shape_fk[dim_idx];
  }

  /* scaler */
  config->scaler =
      AGD_ONE / AGD_SQRT((AGD_TYPE)(config->n_elements_padded_window));

  /* use_taper */
  config->use_taper = (int *)malloc((size_t)(2 * n_dims) * sizeof(int));
  if (config->use_taper == NULL) goto err;
  memcpy(config->use_taper, use_taper,
         (size_t)(2 * (n_dims - 1)) * sizeof(int));
  config->use_taper[2 * (n_dims - 1) + 0] = 1; /* time dimension */
  config->use_taper[2 * (n_dims - 1) + 1] = 1;

  /* taper_windows */
  config->taper_windows =
      (AGD_TYPE **)calloc((size_t)n_dims, sizeof(AGD_TYPE *));
  if (config->taper_windows == NULL) goto err;
  for (dim_idx = 0; dim_idx < n_dims; dim_idx++) {
    if (set_taper_window(config->window_shape[dim_idx],
                         config->use_taper + dim_idx * 2,
                         config->taper_windows + dim_idx))
      goto err;
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_fk_window_config\n");
  return 1;
}

static int append_fk_config(int const n_times, int const n_times_out_capacity,
                            int const n_times_window, int const n_dims,
                            int const *const window_shape,
                            int const *const use_taper,
                            struct FKConfigs *const fk_configs) {
  struct FKWindowConfig *ptr;

  ptr = (struct FKWindowConfig *)realloc(
      fk_configs->window_configs, (size_t)(fk_configs->n_window_configs + 1) *
                                      sizeof(struct FKWindowConfig));
  if (ptr == NULL) goto err;
  fk_configs->window_configs = ptr;
  memset(fk_configs->window_configs + fk_configs->n_window_configs, 0,
         sizeof(struct FKWindowConfig));

  if (set_fk_window_config(
          n_times, n_times_out_capacity, n_times_window, n_dims, window_shape,
          use_taper, fk_configs->window_configs + fk_configs->n_window_configs))
    goto err;

  fk_configs->n_window_configs++;

  return 0;
err:
  fprintf(stderr, "ERROR in append_fk_config\n");
  return 1;
}

static int set_fk_configs(int const n_times, int const n_times_out_capacity,
                          int const n_times_window, int const n_dims,
                          int const *const window_shape,
                          int const *const use_taper,
                          struct FKConfigs *const fk_configs,
                          int *const window_fk_config_idx) {
  int config_idx;

  for (config_idx = 0; config_idx < fk_configs->n_window_configs;
       config_idx++) {
    struct FKWindowConfig const *const config =
        fk_configs->window_configs + config_idx;
    if ((config->n_dims == n_dims) &&
        (config->window_shape[n_dims - 1] == n_times_window) &&
        (get_combined_time(config->n_time_windows,
                           config->window_shape[n_dims - 1]) == n_times)) {
      int dim_idx;
      for (dim_idx = 0; dim_idx < n_dims - 1; dim_idx++) {
        if ((config->window_shape[dim_idx] != window_shape[dim_idx]) ||
            (config->use_taper[dim_idx] != use_taper[dim_idx])) {
          break;
        }
      }
      if (dim_idx == n_dims) {
        *window_fk_config_idx = config_idx;
        return 0;
      }
    }
  }

  /* No existing config matches requirements, so make a new one */
  if (append_fk_config(n_times, n_times_out_capacity, n_times_window, n_dims,
                       window_shape, use_taper, fk_configs))
    goto err;
  *window_fk_config_idx = fk_configs->n_window_configs - 1;

  return 0;
err:
  fprintf(stderr, "ERROR in set_fk_configs\n");
  return 1;
}

static int set_fk_fftw_plans(void *const temp_1, void *const temp_2,
                             struct FKConfigs *const fk_configs) {
  int config_idx;
  for (config_idx = 0; config_idx < fk_configs->n_window_configs;
       config_idx++) {
    struct FKWindowConfig *const config =
        fk_configs->window_configs + config_idx;
    int const rank = config->n_dims;
    int const *const n = config->padded_window_shape;
    int const howmany = config->n_time_windows;
    int const istride = 1;
    int const ostride = 1;
    unsigned flags = AGD_FFTW_FLAGS;
    int const fkdist = config->n_elements_padded_window_fk;
    int const txdist = config->n_elements_padded_window;
    config->fwd_plan = AGD_FFTW_PLAN_C2R(
        rank, n, howmany, (AGD_FFTW_COMPLEX *)temp_1, NULL, istride, fkdist,
        (AGD_TYPE *)temp_2, NULL, ostride, txdist, flags);
    if (config->fwd_plan == NULL) goto err;
    config->adj_plan = AGD_FFTW_PLAN_R2C(
        rank, n, howmany, (AGD_TYPE *)temp_2, NULL, ostride, txdist,
        (AGD_FFTW_COMPLEX *)temp_1, NULL, istride, fkdist, flags);
    if (config->adj_plan == NULL) goto err;
  }

  return 0;
err:
  fprintf(stderr, "ERROR in create_fk_fftw_plans\n");
  return 1;
}

static void free_fk_config(struct FKConfigs *const fk_configs) {
  if (fk_configs->window_configs != NULL) {
    int config_idx;
    for (config_idx = 0; config_idx < fk_configs->n_window_configs;
         config_idx++) {
      struct FKWindowConfig *const config =
          fk_configs->window_configs + config_idx;
      free(config->window_shape);
      config->window_shape = NULL;
      free(config->padded_window_shape);
      config->padded_window_shape = NULL;
      free(config->padded_window_shape_fk);
      config->padded_window_shape_fk = NULL;
      free(config->use_taper);
      config->use_taper = NULL;
      if (config->taper_windows != NULL) {
        int dim_idx;
        for (dim_idx = 0; dim_idx < config->n_dims; dim_idx++) {
          free(config->taper_windows[dim_idx]);
          config->taper_windows[dim_idx] = NULL;
        }
        free(config->taper_windows);
        config->taper_windows = NULL;
      }
      AGD_FFTW_DESTROY_PLAN(config->fwd_plan);
      config->fwd_plan = NULL;
      AGD_FFTW_DESTROY_PLAN(config->adj_plan);
      config->adj_plan = NULL;
    }
    free(fk_configs->window_configs);
    fk_configs->window_configs = NULL;
    fk_configs->n_window_configs = 0;
  }
}

static void fk_forward(AGD_FFTW_COMPLEX const *const fk_data,
                       struct FKWindowConfig const *const config,
                       void *const temp_2, void *const temp_1) {
  int const n_times_window = config->window_shape[config->n_dims - 1];
  int const n_times_padded_window =
      config->padded_window_shape[config->n_dims - 1];
  int const n_times_halfwindow = n_times_window / 2;
  int const n_times_output_capacity = config->n_times_output_capacity;
  AGD_TYPE const *const time_scaler = config->taper_windows[config->n_dims - 1];
  int coords[AGD_MAX_N_COORDS] = ZERO_INIT;
  int trace_idx;

  /* Copy input data into temp_1 */
  memcpy(
      temp_1, fk_data,
      (size_t)(config->n_time_windows * config->n_elements_padded_window_fk) *
          sizeof(AGD_FFTW_COMPLEX));

  /* FK->TX transform - temp_1 -> temp_2 */
  AGD_FFTW_EXECUTE_C2R(config->fwd_plan, (AGD_FFTW_COMPLEX *)temp_1,
                       (AGD_TYPE *)temp_2);

  /* Sum time windows, scale, and remove padding - temp_2 to temp_1 */
  memset(temp_1, 0,
         (size_t)config->n_elements_tx_output_capacity * sizeof(AGD_TYPE));

  for (trace_idx = 0; trace_idx < config->n_traces; trace_idx++) {
    AGD_TYPE scaler = config->scaler;
    int padded_trace_idx = 0;
    int dim_stride = 1;
    int dim_idx;
    int time_window_idx;
    AGD_TYPE *const temp_1_trace =
        (AGD_TYPE *)temp_1 + trace_idx * n_times_output_capacity;
    AGD_TYPE *temp_2_trace;

    for (dim_idx = config->n_dims - 2; dim_idx >= 0; dim_idx--) {
      padded_trace_idx += coords[dim_idx] * dim_stride;
      dim_stride *= config->padded_window_shape[dim_idx];
      scaler *= config->taper_windows[dim_idx][coords[dim_idx]];
    }
    temp_2_trace =
        (AGD_TYPE *)temp_2 + padded_trace_idx * n_times_padded_window;

    for (time_window_idx = 0; time_window_idx < config->n_time_windows;
         time_window_idx++) {
      int time_idx;
      AGD_TYPE *const temp_1_window =
          temp_1_trace + time_window_idx * n_times_halfwindow;
      AGD_TYPE *const temp_2_window =
          temp_2_trace + time_window_idx * config->n_elements_padded_window;
      for (time_idx = 0; time_idx < n_times_window; time_idx++) {
        temp_1_window[time_idx] +=
            temp_2_window[time_idx] * scaler * time_scaler[time_idx];
      }
    }

    for (dim_idx = config->n_dims - 2; dim_idx >= 0; dim_idx--) {
      coords[dim_idx]++;
      if (coords[dim_idx] == config->window_shape[dim_idx]) {
        coords[dim_idx] = 0;
      } else {
        break;
      }
    }
  }
}

static void fk_adjoint(struct FKWindowConfig const *const config,
                       void *const temp_2, void *const temp_1) {
  int const n_times_window = config->window_shape[config->n_dims - 1];
  int const n_times_padded_window =
      config->padded_window_shape[config->n_dims - 1];
  int const n_times_halfwindow = n_times_window / 2;
  int const n_freqs = config->padded_window_shape_fk[config->n_dims - 1];
  int const n_times_output_capacity = config->n_times_output_capacity;
  AGD_TYPE const *const time_scaler = config->taper_windows[config->n_dims - 1];
  int trace_idx;
  int time_window_idx;
  int coords[AGD_MAX_N_COORDS] = ZERO_INIT;

  /* Extract time windows, scale, and add padding, temp_1 -> temp_2 */
  memset(temp_2, 0,
         (size_t)(config->n_time_windows * config->n_elements_padded_window) *
             sizeof(AGD_TYPE));
  for (trace_idx = 0; trace_idx < config->n_traces; trace_idx++) {
    AGD_TYPE scaler = config->scaler;
    int padded_trace_idx = 0;
    int dim_stride = 1;
    int dim_idx;
    AGD_TYPE const *const temp_1_trace =
        (AGD_TYPE *)temp_1 + trace_idx * n_times_output_capacity;
    AGD_TYPE *temp_2_trace;

    for (dim_idx = config->n_dims - 2; dim_idx >= 0; dim_idx--) {
      padded_trace_idx += coords[dim_idx] * dim_stride;
      dim_stride *= config->padded_window_shape[dim_idx];
      scaler *= config->taper_windows[dim_idx][coords[dim_idx]];
    }
    temp_2_trace =
        (AGD_TYPE *)temp_2 + padded_trace_idx * n_times_padded_window;

    for (time_window_idx = 0; time_window_idx < config->n_time_windows;
         time_window_idx++) {
      int time_idx;
      AGD_TYPE const *const temp_1_window =
          temp_1_trace + time_window_idx * n_times_halfwindow;
      AGD_TYPE *const temp_2_window =
          temp_2_trace + time_window_idx * config->n_elements_padded_window;
      for (time_idx = 0; time_idx < n_times_window; time_idx++) {
        temp_2_window[time_idx] =
            temp_1_window[time_idx] * scaler * time_scaler[time_idx];
      }
    }

    for (dim_idx = config->n_dims - 2; dim_idx >= 0; dim_idx--) {
      coords[dim_idx]++;
      if (coords[dim_idx] == config->window_shape[dim_idx]) {
        coords[dim_idx] = 0;
      } else {
        break;
      }
    }
  }

  /* TX->FK transform, temp_2 -> temp_1 */
  AGD_FFTW_EXECUTE_R2C(config->adj_plan, (AGD_TYPE *)temp_2,
                       (AGD_FFTW_COMPLEX *)temp_1);

  /* Multiply elements, except the first and (if number of time samples is
   * even) last in the frequency dimension, by 2 */
  for (time_window_idx = 0; time_window_idx < config->n_time_windows;
       time_window_idx++) {
    AGD_FFTW_COMPLEX *const window =
        (AGD_FFTW_COMPLEX *)temp_1 +
        time_window_idx * config->n_elements_padded_window_fk;
    for (trace_idx = 0; trace_idx < config->n_traces_padded; trace_idx++) {
      int idx;
      AGD_FFTW_COMPLEX *const trace = window + trace_idx * n_freqs;
      for (idx = 2 * 1; idx < 2 * (1 + (n_times_padded_window - 1) / 2);
           idx++) {
        ((AGD_TYPE *)trace)[idx] *= AGD_TWO;
      }
    }
  }
}
