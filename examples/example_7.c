#include <stdio.h>
#include <stdlib.h>

#include "agdeblend.h"

int main(void) {
  size_t n_traces = 64;
  size_t n_times = 1024;
  int n_patches = 1;
  int volumes[1] = {0};
  int n_dims[1] = {2};
  int wavelet_lengths[1] = {513};
  int **window_shapes = (int **)malloc(sizeof(int *));
  int **coords = (int **)malloc(sizeof(int *));
  int **shapes = (int **)malloc(sizeof(int *));
  int long **shottimes = (int long **)malloc(sizeof(int long *));
  int **channels = (int **)malloc(sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(sizeof(enum AGDTraceType *));
  float **data = (float **)malloc(sizeof(float *));
  float **wavelets = (float **)malloc(sizeof(float *));
  int **wavelet_idxs = (int **)malloc(sizeof(int *));
  float initial_factor = 1.0f;
  int n_its = 1250;
  int print_freq = -1;
  FILE *fid;

  window_shapes[0] = (int *)malloc(n_dims[0] * sizeof(int));
  coords[0] = (int *)malloc((n_dims[0] - 1) * sizeof(int));
  shapes[0] = (int *)malloc(n_dims[0] * sizeof(int));
  shottimes[0] = (int long *)malloc(n_traces * sizeof(int long));
  channels[0] = (int *)malloc(n_traces * sizeof(int));
  trace_types[0] =
      (enum AGDTraceType *)malloc(n_traces * sizeof(enum AGDTraceType));
  data[0] = (float *)malloc(n_traces * n_times * sizeof(float));
  wavelets[0] = (float *)malloc(wavelet_lengths[0] * sizeof(float));
  wavelet_idxs[0] = (int *)malloc(n_traces * sizeof(int));

  window_shapes[0][0] = 32;
  window_shapes[0][1] = 256;
  coords[0][0] = 0;
  shapes[0][0] = n_traces;
  shapes[0][1] = n_times;

  /* Load from input files */
  fid = fopen("data_2_shottimes.bin", "rb");
  fread(shottimes[0], sizeof(long int), n_traces, fid);
  fclose(fid);

  fid = fopen("data_2_channels.bin", "rb");
  fread(channels[0], sizeof(int), n_traces, fid);
  fclose(fid);

  fid = fopen("data_2_trace_types.bin", "rb");
  fread(trace_types[0], sizeof(enum AGDTraceType), n_traces, fid);
  fclose(fid);

  fid = fopen("data_2_blended_data.bin", "rb");
  fread(data[0], sizeof(float), n_traces * n_times, fid);
  fclose(fid);

  fid = fopen("data_2_wavelet.bin", "rb");
  fread(wavelets[0], sizeof(float), wavelet_lengths[0], fid);
  fclose(fid);

  fid = fopen("data_2_wavelet_idxs.bin", "rb");
  fread(wavelet_idxs[0], sizeof(int), n_traces, fid);
  fclose(fid);

  /* Deblend */
  agd_deblend(n_patches, volumes, n_dims, (int const *const *)window_shapes,
              (int const *const *)coords, (int const *const *)shapes,
              (long int const *const *)shottimes, (int const *const *)channels,
              (enum AGDTraceType const *const *)trace_types, wavelet_lengths,
              (int const *const *)wavelet_idxs, (float const *const *)wavelets,
              initial_factor, n_its, print_freq, data);

  /* Write to file */
  fid = fopen("out/data_2_deblended_data.bin", "wb");
  fwrite(data[0], sizeof(float), n_traces * n_times, fid);
  fclose(fid);

  free(window_shapes[0]);
  free(coords[0]);
  free(shapes[0]);
  free(shottimes[0]);
  free(channels[0]);
  free(trace_types[0]);
  free(data[0]);
  free(wavelets[0]);
  free(wavelet_idxs[0]);
  free(window_shapes);
  free(coords);
  free(shapes);
  free(shottimes);
  free(channels);
  free(trace_types);
  free(data);
  free(wavelets);
  free(wavelet_idxs);

  return 0;
}
