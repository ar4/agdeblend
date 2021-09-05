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

static int set_wavelet_idxs(int const n_traces, int const n_times,
                            enum AGDTraceType const *const trace_types,
                            int const *const wavelet_idxs,
                            struct WaveletConfig *const config,
                            struct WaveletParams *const params) {
  int trace_idx;
  int const n_freqs = n_times / 2 + 1;

  params->source_config_idxs = (int *)malloc((size_t)n_traces * sizeof(int));
  if (params->source_config_idxs == NULL) goto err;

  for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
    int const wavelet_idx = wavelet_idxs[trace_idx];
    if (trace_types[trace_idx] == AGDLive) {
      int source_config_idx;
      for (source_config_idx = 0; source_config_idx < config->n_source_configs;
           source_config_idx++) {
        struct WaveletSourceConfig const *const source_config =
            config->source_configs + source_config_idx;
        if (source_config->original_idx == wavelet_idx &&
            source_config->n_times == n_times) {
          /* A suitable source config already exists */
          params->source_config_idxs[trace_idx] = source_config_idx;
          break;
        }
      }
      if (source_config_idx == config->n_source_configs) {
        /* Create a new source config */
        struct WaveletSourceConfig *const ptr =
            (struct WaveletSourceConfig *)realloc(
                config->source_configs, (size_t)(config->n_source_configs + 1) *
                                            sizeof(struct WaveletSourceConfig));
        if (ptr == NULL) goto err;
        config->source_configs = ptr;
        config->n_source_configs++;
        params->source_config_idxs[trace_idx] = source_config_idx;
        config->source_configs[source_config_idx].original_idx = wavelet_idx;
        config->source_configs[source_config_idx].n_times = n_times;
        config->source_configs[source_config_idx].wavelet =
            (AGD_FFTW_COMPLEX *)AGD_FFTW_ALLOC_COMPLEX((size_t)n_freqs);
      }
    } else {
      /* Missing or bad trace */
      params->source_config_idxs[trace_idx] = -1;
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_wavelet_idxs\n");
  return 1;
}

static int set_wavelet_plan(int const n_traces, int const n_times,
                            struct WaveletConfig *const config,
                            struct WaveletParams *const params) {
  int plan_idx;

  for (plan_idx = 0; plan_idx < config->n_fftw_plans; plan_idx++) {
    struct WaveletFFTWPlan const *const plan = config->plans + plan_idx;
    if (plan->n_traces == n_traces && plan->n_times == n_times) {
      /* Suitable FFTW plan already exists */
      params->plan_idx = plan_idx;
      break;
    }
  }

  if (plan_idx == config->n_fftw_plans) {
    /* Add new FFTW plan */
    struct WaveletFFTWPlan *const ptr = (struct WaveletFFTWPlan *)realloc(
        config->plans,
        (size_t)(config->n_fftw_plans + 1) * sizeof(struct WaveletFFTWPlan));
    if (ptr == NULL) goto err;
    config->plans = ptr;
    config->n_fftw_plans++;
    params->plan_idx = plan_idx;
    config->plans[plan_idx].n_traces = n_traces;
    config->plans[plan_idx].n_times = n_times;
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_wavelet_plan\n");
  return 1;
}

static int set_wavelet_params(int const n_traces, int const n_times,
                              enum AGDTraceType const *const trace_types,
                              int const *const wavelet_idxs,
                              struct WaveletConfig *const config,
                              struct WaveletParams *const params) {
  if (wavelet_idxs == NULL) return 0;

  if (set_wavelet_idxs(n_traces, n_times, trace_types, wavelet_idxs, config,
                       params))
    goto err;

  if (set_wavelet_plan(n_traces, n_times, config, params)) goto err;

  return 0;
err:
  fprintf(stderr, "ERROR in set_wavelet_params\n");
  return 1;
}

static void free_wavelet_config(struct WaveletConfig *const config) {
  if (config->source_configs != NULL) {
    int idx;
    for (idx = 0; idx < config->n_source_configs; idx++) {
      struct WaveletSourceConfig *const source_config =
          config->source_configs + idx;
      AGD_FFTW_FREE(source_config->wavelet);
      source_config->wavelet = NULL;
    }
    free(config->source_configs);
    config->source_configs = NULL;
  }

  if (config->plans != NULL) {
    int idx;
    for (idx = 0; idx < config->n_fftw_plans; idx++) {
      struct WaveletFFTWPlan *const plan = config->plans + idx;
      AGD_FFTW_DESTROY_PLAN(plan->forward_plan);
      AGD_FFTW_DESTROY_PLAN(plan->backward_plan);
      plan->forward_plan = NULL;
      plan->backward_plan = NULL;
    }
    free(config->plans);
    config->plans = NULL;
  }
  config->n_source_configs = 0;
  config->n_fftw_plans = 0;
}

static void free_wavelet_params(struct WaveletParams *const params) {
  free(params->source_config_idxs);
  params->source_config_idxs = NULL;
}

static int set_wavelet_fftw_plans(void *const temp_1, void *const temp_2,
                                  struct WaveletConfig *const config) {
  int plan_idx;
  for (plan_idx = 0; plan_idx < config->n_fftw_plans; plan_idx++) {
    struct WaveletFFTWPlan *const plan = config->plans + plan_idx;
    int const rank = 1;
    int n[1];
    int const howmany = plan->n_traces;
    int const istride = 1;
    int const ostride = 1;
    unsigned flags = AGD_FFTW_FLAGS;
    int const fkdist = plan->n_times / 2 + 1;
    int const txdist = plan->n_times;
    n[0] = plan->n_times;
    plan->forward_plan = AGD_FFTW_PLAN_R2C(
        rank, n, howmany, (AGD_TYPE *)temp_1, NULL, istride, txdist,
        (AGD_FFTW_COMPLEX *)temp_2, NULL, ostride, fkdist, flags);
    if (plan->forward_plan == NULL) goto err;
    plan->backward_plan = AGD_FFTW_PLAN_C2R(
        rank, n, howmany, (AGD_FFTW_COMPLEX *)temp_2, NULL, ostride, fkdist,
        (AGD_TYPE *)temp_1, NULL, istride, txdist, flags);
    if (plan->backward_plan == NULL) goto err;
  }

  return 0;
err:
  fprintf(stderr, "ERROR in create_wavelet_fftw_plans\n");
  return 1;
}

static int transform_wavelets(struct WaveletConfig *const config,
                              int const *const wavelet_lengths,
                              AGD_TYPE const *const *const wavelets,
                              void *const temp_1) {
  int source_config_idx;
  for (source_config_idx = 0; source_config_idx < config->n_source_configs;
       source_config_idx++) {
    struct WaveletSourceConfig *const source_config =
        config->source_configs + source_config_idx;
    AGD_TYPE const *const original_wavelet =
        wavelets[source_config->original_idx];
    int const rank = 1;
    int n[1];
    int const howmany = 1;
    int const istride = 1;
    int const ostride = 1;
    unsigned const flags = FFTW_ESTIMATE;
    int const fkdist = source_config->n_times / 2 + 1;
    int const txdist = source_config->n_times;
    AGD_TYPE const scale = AGD_ONE / (AGD_TYPE)source_config->n_times;
    AGD_FFTW_PLAN plan;
    int idx;
    n[0] = txdist;
    plan = AGD_FFTW_PLAN_R2C(rank, n, howmany, (AGD_TYPE *)temp_1, NULL,
                             istride, txdist, source_config->wavelet, NULL,
                             ostride, fkdist, flags);
    if (plan == NULL) goto err;

    memset(temp_1, 0, (size_t)txdist * sizeof(AGD_TYPE));
    memcpy(temp_1, original_wavelet,
           (size_t)wavelet_lengths[source_config->original_idx] *
               sizeof(AGD_TYPE));
    AGD_FFTW_EXECUTE(plan);
    AGD_FFTW_DESTROY_PLAN(plan);

    /* Convert to have norm equal to 1/n_times */
    for (idx = 0; idx < fkdist; idx++) {
      AGD_TYPE const re = source_config->wavelet[idx][0];
      AGD_TYPE const im = source_config->wavelet[idx][1];
      AGD_TYPE const angle = AGD_ATAN2(im, re);
      source_config->wavelet[idx][0] = AGD_COS(angle) * scale;
      source_config->wavelet[idx][1] = AGD_SIN(angle) * scale;
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in transform_wavelets\n");
  return 1;
}

static void wavelet_operator(struct WaveletConfig const *const config,
                             struct WaveletParams const *const params,
                             void *const temp_1, void *const temp_2,
                             AGD_TYPE const direction) {
  struct WaveletFFTWPlan const *const plan = config->plans + params->plan_idx;
  int n_times = -1;
  int n_freqs;
  int trace_idx;

  if (config->n_fftw_plans == 0) return; /* Not using wavelets */

  /* Find n_times from a non-missing trace (if one exists) */
  for (trace_idx = 0; trace_idx < plan->n_traces; trace_idx++) {
    int const source_config_idx = params->source_config_idxs[trace_idx];
    if (source_config_idx >= 0) {
      n_times = config->source_configs[source_config_idx].n_times;
      break;
    }
  }
  if (n_times < 0) return; /* No non-missing traces */
  n_freqs = n_times / 2 + 1;

  /* FFT input, temp_1 -> temp_2 */
  AGD_FFTW_EXECUTE_R2C(plan->forward_plan, (AGD_TYPE *)temp_1,
                       (AGD_FFTW_COMPLEX *)temp_2);

  /* Multiply by wavelet */
  for (trace_idx = 0; trace_idx < plan->n_traces; trace_idx++) {
    int const source_config_idx = params->source_config_idxs[trace_idx];
    struct WaveletSourceConfig const *const source_config =
        config->source_configs + source_config_idx;
    AGD_FFTW_COMPLEX const *wavelet;
    AGD_FFTW_COMPLEX *trace;
    int freq_idx;
    if (source_config_idx < 0) continue;
    wavelet = (AGD_FFTW_COMPLEX const *)source_config->wavelet;
    trace = (AGD_FFTW_COMPLEX *)temp_2 + trace_idx * n_freqs;
    for (freq_idx = 0; freq_idx < n_freqs; freq_idx++) {
      AGD_TYPE const re = trace[freq_idx][0];
      AGD_TYPE const im = trace[freq_idx][1];
      AGD_TYPE const wre = wavelet[freq_idx][0];
      AGD_TYPE const wim = wavelet[freq_idx][1];
      trace[freq_idx][0] = re * wre - direction * im * wim;
      trace[freq_idx][1] = direction * re * wim + im * wre;
    }
  }

  /* Inverse FFT, temp_2 -> temp_1 */
  AGD_FFTW_EXECUTE_C2R(plan->backward_plan, (AGD_FFTW_COMPLEX *)temp_2,
                       (AGD_TYPE *)temp_1);
}

static void wavelet_forward(struct WaveletConfig const *const config,
                            struct WaveletParams const *const params,
                            void *const temp_1, void *const temp_2) {
  wavelet_operator(config, params, temp_1, temp_2, +AGD_ONE);
}

static void wavelet_adjoint(struct WaveletConfig const *const config,
                            struct WaveletParams const *const params,
                            void *const temp_1, void *const temp_2) {
  wavelet_operator(config, params, temp_1, temp_2, -AGD_ONE);
}
