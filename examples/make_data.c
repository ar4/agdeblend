#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "agdeblend.h"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

/* 1 2D volume */
static void make_1(void) {
  size_t n_traces = 128;
  size_t n_times = 512;
  size_t n_times_out = 768;
  int long delay = 256;
  int long jitter = 50;
  int n_traces_arr[1];
  int n_times_arr[1];
  int n_times_out_arr[1];
  int long **shottimes = (int long **)malloc(sizeof(int long *));
  int long **shottimes_out = (int long **)malloc(sizeof(int long *));
  int **channels = (int **)malloc(sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(sizeof(enum AGDTraceType *));
  float true_trace[512];
  fftwf_complex true_trace_f[512 / 2 + 1];
  float **true_data = (float **)malloc(sizeof(float *));
  float **blended_data = (float **)malloc(sizeof(float *));
  float **blended_data2 = (float **)malloc(sizeof(float *));
  size_t trace_idx;
  size_t time_idx;
  fftwf_plan r2c_plan;
  fftwf_plan c2r_plan;
  FILE *fid;

  n_traces_arr[0] = (int)n_traces;
  n_times_arr[0] = (int)n_times;
  n_times_out_arr[0] = (int)n_times_out;

  shottimes[0] = (int long *)malloc(n_traces * sizeof(int long));
  shottimes_out[0] = (int long *)malloc(n_traces * sizeof(int long));
  channels[0] = (int *)calloc(n_traces, sizeof(int));
  trace_types[0] =
      (enum AGDTraceType *)calloc(n_traces, sizeof(enum AGDTraceType));
  true_data[0] = (float *)malloc(n_traces * n_times * sizeof(float));
  blended_data[0] = (float *)malloc(n_traces * n_times * sizeof(float));
  blended_data2[0] = (float *)malloc(n_traces * n_times_out * sizeof(float));

  srand(1);

  /* Set shottimes */
  for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
    shottimes[0][trace_idx] =
        (long int)trace_idx * delay +
        (long int)((((float)rand()) / ((float)RAND_MAX) - 0.5f) * 2.0f *
                   (float)jitter);
  }

  /* Set shottimes_out to be the same but earlier */
  for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
    shottimes_out[0][trace_idx] = shottimes[0][trace_idx] - 128;
  }

  /* Set data */
  /* Begin by creating a random trace and high-pass filtering it
   * so that the maximum wavelength is 128 samples. This is then
   * used for all traces. */
  for (time_idx = 0; time_idx < n_times; time_idx++) {
    true_trace[time_idx] = (((float)rand()) / ((float)RAND_MAX) - 0.5f) * 2.0f;
  }
  r2c_plan = fftwf_plan_dft_r2c_1d((int)n_times, true_trace, true_trace_f,
                                   FFTW_ESTIMATE);
  c2r_plan = fftwf_plan_dft_c2r_1d((int)n_times, true_trace_f, true_trace,
                                   FFTW_ESTIMATE);
  fftwf_execute(r2c_plan);
  memset(true_trace_f, 0, 4 * sizeof(fftwf_complex));
  fftwf_execute(c2r_plan);
  fftwf_destroy_plan(r2c_plan);
  fftwf_destroy_plan(c2r_plan);

  for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
    memcpy(true_data[0] + trace_idx * n_times, true_trace,
           n_times * sizeof(float));
  }

  /* Blend */
  agd_blend(1, n_traces_arr, n_times_arr, (long int const *const *)shottimes,
            (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, true_data,
            AGDBlendSum, 0, 1, n_traces_arr, n_times_arr,
            (long int const *const *)shottimes, (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, blended_data);

  /* Blend and output with earlier shottimes and longer traces */
  agd_blend(1, n_traces_arr, n_times_arr, (long int const *const *)shottimes,
            (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, true_data,
            AGDBlendSum, 0, 1, n_traces_arr, n_times_out_arr,
            (long int const *const *)shottimes_out,
            (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, blended_data2);

  /* Write to files */
  fid = fopen("data_1_shottimes.bin", "wb");
  fwrite(shottimes[0], sizeof(long int), n_traces, fid);
  fclose(fid);

  fid = fopen("data_1_shottimes_out.bin", "wb");
  fwrite(shottimes_out[0], sizeof(long int), n_traces, fid);
  fclose(fid);

  fid = fopen("data_1_channels.bin", "wb");
  fwrite(channels[0], sizeof(int), n_traces, fid);
  fclose(fid);

  fid = fopen("data_1_trace_types.bin", "wb");
  fwrite(trace_types[0], sizeof(enum AGDTraceType), n_traces, fid);
  fclose(fid);

  fid = fopen("data_1_true_data.bin", "wb");
  fwrite(true_data[0], sizeof(float), n_traces * n_times, fid);
  fclose(fid);

  fid = fopen("data_1_blended_data.bin", "wb");
  fwrite(blended_data[0], sizeof(float), n_traces * n_times, fid);
  fclose(fid);

  fid = fopen("data_1_blended_data_2.bin", "wb");
  fwrite(blended_data2[0], sizeof(float), n_traces * n_times_out, fid);
  fclose(fid);

  free(shottimes[0]);
  free(shottimes_out[0]);
  free(channels[0]);
  free(trace_types[0]);
  free(true_data[0]);
  free(blended_data[0]);
  free(blended_data2[0]);
  free(shottimes);
  free(shottimes_out);
  free(channels);
  free(trace_types);
  free(true_data);
  free(blended_data);
  free(blended_data2);
}

/* 2D Volume with chirp wavelet */
static void make_2(void) {
  size_t n_traces = 64;
  size_t n_times_preconv = 512;
  size_t n_times_wavelet = 513;
  size_t n_times = n_times_preconv + n_times_wavelet - 1;
  size_t n_freqs = n_times / 2 + 1;
  int long delay = 256;
  int long jitter = 50;
  float chirp_start_freq = 0.01f;
  float chirp_end_freq = 0.15f;
  float beta = (chirp_end_freq - chirp_start_freq) / n_times_wavelet;
  int n_chirp_taper = (int)(0.5f / chirp_start_freq);
  int n_traces_arr[1];
  int n_times_arr[1];
  int long **shottimes = (int long **)malloc(sizeof(int long *));
  int **channels = (int **)malloc(sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(sizeof(enum AGDTraceType *));
  float **true_data = (float **)malloc(sizeof(float *));
  float **blended_data = (float **)malloc(sizeof(float *));
  float *wavelet = (float *)calloc(n_times, sizeof(float));
  int **wavelet_idxs = (int **)malloc(sizeof(int *));
  float *true_trace = (float *)calloc(n_times, sizeof(float));
  fftwf_complex *true_trace_f =
      (fftwf_complex *)fftwf_malloc(n_freqs * sizeof(fftwf_complex));
  fftwf_complex *wavelet_f =
      (fftwf_complex *)fftwf_malloc(n_freqs * sizeof(fftwf_complex));
  size_t trace_idx;
  size_t time_idx;
  fftwf_plan r2c_plan;
  fftwf_plan c2r_plan;
  FILE *fid;

  n_traces_arr[0] = (int)n_traces;
  n_times_arr[0] = (int)(n_times);

  shottimes[0] = (int long *)malloc(n_traces * sizeof(int long));
  channels[0] = (int *)calloc(n_traces, sizeof(int));
  trace_types[0] =
      (enum AGDTraceType *)calloc(n_traces, sizeof(enum AGDTraceType));
  true_data[0] = (float *)calloc(n_traces * n_times, sizeof(float));
  blended_data[0] = (float *)malloc(n_traces * n_times * sizeof(float));
  wavelet_idxs[0] = (int *)calloc(n_traces, sizeof(int));

  srand(1);

  /* Set shottimes */
  for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
    shottimes[0][trace_idx] =
        (long int)trace_idx * delay +
        (long int)((((float)rand()) / ((float)RAND_MAX) - 0.5f) * 2.0f *
                   (float)jitter);
  }

  /* Set wavelet */
  for (time_idx = 0; time_idx < n_times_wavelet; time_idx++) {
    float phase =
        2.0f * M_PI *
        (chirp_start_freq * time_idx + 0.5f * beta * time_idx * time_idx);
    int dist_from_edge =
        time_idx < n_times_wavelet / 2 ? time_idx : n_times_wavelet - time_idx;
    float taper = dist_from_edge < n_chirp_taper
                      ? sinf(M_PI / 2.0f * dist_from_edge / n_chirp_taper)
                      : 1.0f;
    wavelet[time_idx] = cosf(phase) * taper;
  }

  /* Set data */
  /* Begin by creating a random trace, convolving it with the source wavelet,
   * and high-pass filtering it so that the maximum wavelength is 128 samples.
   * This is then used for all traces. */
  for (time_idx = 0; time_idx < n_times_preconv; time_idx++) {
    true_trace[time_idx] = (((float)rand()) / ((float)RAND_MAX) - 0.5f) * 2.0f;
  }
  r2c_plan = fftwf_plan_dft_r2c_1d((int)n_times, true_trace, true_trace_f,
                                   FFTW_ESTIMATE);
  c2r_plan = fftwf_plan_dft_c2r_1d((int)n_times, true_trace_f, true_trace,
                                   FFTW_ESTIMATE);
  fftwf_execute(r2c_plan);
  fftwf_destroy_plan(r2c_plan);
  r2c_plan =
      fftwf_plan_dft_r2c_1d((int)n_times, wavelet, wavelet_f, FFTW_ESTIMATE);
  fftwf_execute(r2c_plan);
  fftwf_destroy_plan(r2c_plan);

  for (time_idx = 0; time_idx < n_freqs; time_idx++) {
    float trace_re = true_trace_f[time_idx][0];
    float trace_im = true_trace_f[time_idx][1];
    float wavelet_re = wavelet_f[time_idx][0];
    float wavelet_im = wavelet_f[time_idx][1];
    true_trace_f[time_idx][0] = trace_re * wavelet_re - trace_im * wavelet_im;
    true_trace_f[time_idx][1] = trace_re * wavelet_im + trace_im * wavelet_re;
  }

  memset(true_trace_f, 0, 4 * sizeof(fftwf_complex));
  fftwf_execute(c2r_plan);
  fftwf_destroy_plan(c2r_plan);

  for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
    memcpy(true_data[0] + trace_idx * n_times, true_trace,
           n_times * sizeof(float));
  }

  /* Blend */
  agd_blend(1, n_traces_arr, n_times_arr, (long int const *const *)shottimes,
            (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, true_data,
            AGDBlendSum, 0, 1, n_traces_arr, n_times_arr,
            (long int const *const *)shottimes, (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, blended_data);

  /* Write to files */
  fid = fopen("data_2_shottimes.bin", "wb");
  fwrite(shottimes[0], sizeof(long int), n_traces, fid);
  fclose(fid);

  fid = fopen("data_2_channels.bin", "wb");
  fwrite(channels[0], sizeof(int), n_traces, fid);
  fclose(fid);

  fid = fopen("data_2_trace_types.bin", "wb");
  fwrite(trace_types[0], sizeof(enum AGDTraceType), n_traces, fid);
  fclose(fid);

  fid = fopen("data_2_true_data.bin", "wb");
  fwrite(true_data[0], sizeof(float), n_traces * n_times, fid);
  fclose(fid);

  fid = fopen("data_2_blended_data.bin", "wb");
  fwrite(blended_data[0], sizeof(float), n_traces * n_times, fid);
  fclose(fid);

  fid = fopen("data_2_wavelet.bin", "wb");
  fwrite(wavelet, sizeof(float), n_times_wavelet, fid);
  fclose(fid);

  fid = fopen("data_2_wavelet_idxs.bin", "wb");
  fwrite(wavelet_idxs[0], sizeof(int), n_traces, fid);
  fclose(fid);

  free(shottimes[0]);
  free(channels[0]);
  free(trace_types[0]);
  free(true_data[0]);
  free(blended_data[0]);
  free(wavelet_idxs[0]);
  free(shottimes);
  free(channels);
  free(trace_types);
  free(true_data);
  free(blended_data);
  free(wavelet);
  free(wavelet_idxs);
  free(true_trace);
  fftwf_free(true_trace_f);
  fftwf_free(wavelet_f);
}

/* 2 2D volumes */
static void make_3() {
  size_t n_traces = 128;
  size_t n_times[2] = {512, 768};
  int long delay = 256;
  int long jitter = 50;
  int n_traces_arr[2];
  int n_times_arr[2];
  int long **shottimes = (int long **)malloc(2 * sizeof(int long *));
  int **channels = (int **)malloc(2 * sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(2 * sizeof(enum AGDTraceType *));
  float true_trace[768];
  fftwf_complex true_trace_f[768 / 2 + 1];
  float **true_data = (float **)malloc(2 * sizeof(float *));
  float **blended_data = (float **)malloc(2 * sizeof(float *));
  size_t trace_idx;
  size_t time_idx;
  int volume_idx;
  fftwf_plan r2c_plan;
  fftwf_plan c2r_plan;
  char filename[256];
  FILE *fid;

  srand(1);

  for (volume_idx = 0; volume_idx < 2; volume_idx++) {
    n_traces_arr[volume_idx] = (int)n_traces;
    n_times_arr[volume_idx] = (int)n_times[volume_idx];

    shottimes[volume_idx] = (int long *)malloc(n_traces * sizeof(int long));
    channels[volume_idx] = (int *)calloc(n_traces, sizeof(int));
    trace_types[volume_idx] =
        (enum AGDTraceType *)calloc(n_traces, sizeof(enum AGDTraceType));
    true_data[volume_idx] =
        (float *)malloc(n_traces * n_times[volume_idx] * sizeof(float));
    blended_data[volume_idx] =
        (float *)malloc(n_traces * n_times[volume_idx] * sizeof(float));

    /* Set shottimes */
    for (trace_idx = volume_idx; trace_idx < n_traces; trace_idx++) {
      shottimes[volume_idx][trace_idx] =
          (long int)trace_idx * delay +
          (long int)((((float)rand()) / ((float)RAND_MAX) - 0.5f) * 2.0f *
                     (float)jitter);
    }

    /* Set data */
    /* Begin by creating a random trace and high-pass filtering it
     * so that the maximum wavelength is 128 samples. This is then
     * used for all traces. */
    for (time_idx = 0; time_idx < n_times[volume_idx]; time_idx++) {
      true_trace[time_idx] =
          (((float)rand()) / ((float)RAND_MAX) - 0.5f) * 2.0f;
    }
    r2c_plan = fftwf_plan_dft_r2c_1d((int)n_times[volume_idx], true_trace,
                                     true_trace_f, FFTW_ESTIMATE);
    c2r_plan = fftwf_plan_dft_c2r_1d((int)n_times[volume_idx], true_trace_f,
                                     true_trace, FFTW_ESTIMATE);
    fftwf_execute(r2c_plan);
    memset(true_trace_f, 0, 4 * sizeof(fftwf_complex));
    fftwf_execute(c2r_plan);
    fftwf_destroy_plan(r2c_plan);
    fftwf_destroy_plan(c2r_plan);

    for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
      memcpy(true_data[volume_idx] + trace_idx * n_times[volume_idx],
             true_trace, n_times[volume_idx] * sizeof(float));
    }
  }

  /* Blend */
  agd_blend(2, n_traces_arr, n_times_arr, (long int const *const *)shottimes,
            (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, true_data,
            AGDBlendSum, 0, 2, n_traces_arr, n_times_arr,
            (long int const *const *)shottimes, (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, blended_data);

  /* Write to files */
  for (volume_idx = 0; volume_idx < 2; volume_idx++) {
    sprintf(filename, "data_3_shottimes_%d.bin", volume_idx);
    fid = fopen(filename, "wb");
    fwrite(shottimes[volume_idx], sizeof(long int), n_traces, fid);
    fclose(fid);

    sprintf(filename, "data_3_channels_%d.bin", volume_idx);
    fid = fopen(filename, "wb");
    fwrite(channels[volume_idx], sizeof(int), n_traces, fid);
    fclose(fid);

    sprintf(filename, "data_3_trace_types_%d.bin", volume_idx);
    fid = fopen(filename, "wb");
    fwrite(trace_types[volume_idx], sizeof(enum AGDTraceType), n_traces, fid);
    fclose(fid);

    sprintf(filename, "data_3_true_data_%d.bin", volume_idx);
    fid = fopen(filename, "wb");
    fwrite(true_data[volume_idx], sizeof(float), n_traces * n_times[volume_idx],
           fid);
    fclose(fid);

    sprintf(filename, "data_3_blended_data_%d.bin", volume_idx);
    fid = fopen(filename, "wb");
    fwrite(blended_data[volume_idx], sizeof(float),
           n_traces * n_times[volume_idx], fid);
    fclose(fid);

    free(shottimes[volume_idx]);
    free(channels[volume_idx]);
    free(trace_types[volume_idx]);
    free(true_data[volume_idx]);
    free(blended_data[volume_idx]);
  }
  free(shottimes);
  free(channels);
  free(trace_types);
  free(true_data);
  free(blended_data);
}

static void shift_trace(fftwf_complex *trace_f, int n_freqs, float shift,
                        fftwf_complex *trace_f_shifted) {
  int freq_idx;
  for (freq_idx = 0; freq_idx < n_freqs; freq_idx++) {
    float k = (float)freq_idx / (float)(n_freqs - 1) / 2.0f;
    float arg = -2.0f * M_PI * shift * k;
    float shift_re = cosf(arg);
    float shift_im = sinf(arg);
    float trace_re = trace_f[freq_idx][0];
    float trace_im = trace_f[freq_idx][1];
    trace_f_shifted[freq_idx][0] = trace_re * shift_re - trace_im * shift_im;
    trace_f_shifted[freq_idx][1] = trace_re * shift_im + trace_im * shift_re;
  }
}

/* 3D volume, irregular shape */
static void make_4() {
  size_t nx = 16;
  size_t ny = 24;
  size_t n_traces = nx * ny;
  size_t n_times = 512;
  size_t n_freqs = 512 / 2 + 1;
  int long delay = 256;
  int long jitter = 50;
  int n_traces_arr[2];
  int n_times_arr[2];
  int long **shottimes = (int long **)malloc(2 * sizeof(int long *));
  int **channels = (int **)malloc(2 * sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(2 * sizeof(enum AGDTraceType *));
  float true_trace[512];
  fftwf_complex true_trace_f[512 / 2 + 1];
  fftwf_complex true_trace_f_shifted[512 / 2 + 1];
  float **true_data = (float **)malloc(2 * sizeof(float *));
  float **blended_data = (float **)malloc(2 * sizeof(float *));
  size_t patch_idx;
  size_t trace_idx;
  size_t time_idx;
  size_t x;
  size_t y;
  fftwf_plan r2c_plan;
  fftwf_plan c2r_plan;
  char filename[256];
  FILE *fid;

  srand(1);

  for (patch_idx = 0; patch_idx < 2; patch_idx++) {
    n_traces_arr[patch_idx] = (int)n_traces;
    n_times_arr[patch_idx] = (int)n_times;

    shottimes[patch_idx] = (int long *)malloc(n_traces * sizeof(int long));
    channels[patch_idx] = (int *)calloc(n_traces, sizeof(int));
    trace_types[patch_idx] =
        (enum AGDTraceType *)calloc(n_traces, sizeof(enum AGDTraceType));
    true_data[patch_idx] = (float *)malloc(n_traces * n_times * sizeof(float));
    blended_data[patch_idx] =
        (float *)malloc(n_traces * n_times * sizeof(float));

    /* Set shottimes */
    for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
      shottimes[patch_idx][trace_idx] =
          (long int)(patch_idx * n_traces + trace_idx) * delay +
          (long int)((((float)rand()) / ((float)RAND_MAX) - 0.5f) * 2.0f *
                     (float)jitter);
    }
  }

  /* Set data */
  /* Begin by creating a random trace and high-pass filtering it
   * so that the maximum wavelength is 128 samples. This is then
   * used for all traces with a shift that increases with x and y. */
  for (time_idx = 0; time_idx < n_times; time_idx++) {
    true_trace[time_idx] = (((float)rand()) / ((float)RAND_MAX) - 0.5f) * 2.0f;
  }
  r2c_plan = fftwf_plan_dft_r2c_1d((int)n_times, true_trace, true_trace_f,
                                   FFTW_ESTIMATE);
  c2r_plan = fftwf_plan_dft_c2r_1d((int)n_times, true_trace_f_shifted,
                                   true_trace, FFTW_ESTIMATE);
  fftwf_execute(r2c_plan);
  fftwf_destroy_plan(r2c_plan);
  memset(true_trace_f, 0, 4 * sizeof(fftwf_complex));

  for (patch_idx = 0; patch_idx < 2; patch_idx++) {
    for (x = 0; x < nx; x++) {
      for (y = 0; y < ny; y++) {
        shift_trace(true_trace_f, n_freqs,
                    (float)(x + y + patch_idx * (nx + 16)) * 0.2f,
                    true_trace_f_shifted);
        fftwf_execute(c2r_plan);
        memcpy(true_data[patch_idx] + x * ny * n_times + y * n_times,
               true_trace, n_times * sizeof(float));
      }
    }
  }
  fftwf_destroy_plan(c2r_plan);

  /* Blend */
  agd_blend(2, n_traces_arr, n_times_arr, (long int const *const *)shottimes,
            (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, true_data,
            AGDBlendSum, 0, 2, n_traces_arr, n_times_arr,
            (long int const *const *)shottimes, (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, blended_data);

  for (patch_idx = 0; patch_idx < 2; patch_idx++) {
    /* Write to files */
    sprintf(filename, "data_4_shottimes_%ld.bin", patch_idx);
    fid = fopen(filename, "wb");
    fwrite(shottimes[patch_idx], sizeof(long int), n_traces, fid);
    fclose(fid);

    sprintf(filename, "data_4_channels_%ld.bin", patch_idx);
    fid = fopen(filename, "wb");
    fwrite(channels[patch_idx], sizeof(int), n_traces, fid);
    fclose(fid);

    sprintf(filename, "data_4_trace_types_%ld.bin", patch_idx);
    fid = fopen(filename, "wb");
    fwrite(trace_types[patch_idx], sizeof(enum AGDTraceType), n_traces, fid);
    fclose(fid);

    sprintf(filename, "data_4_true_data_%ld.bin", patch_idx);
    fid = fopen(filename, "wb");
    fwrite(true_data[patch_idx], sizeof(float), n_traces * n_times, fid);
    fclose(fid);

    sprintf(filename, "data_4_blended_data_%ld.bin", patch_idx);
    fid = fopen(filename, "wb");
    fwrite(blended_data[patch_idx], sizeof(float), n_traces * n_times, fid);
    fclose(fid);

    free(shottimes[patch_idx]);
    free(channels[patch_idx]);
    free(trace_types[patch_idx]);
    free(true_data[patch_idx]);
    free(blended_data[patch_idx]);
  }
  free(shottimes);
  free(channels);
  free(trace_types);
  free(true_data);
  free(blended_data);
}

int main(void) {
  make_1();
  make_2();
  make_3();
  make_4();
  return 0;
}
