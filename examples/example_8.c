#include <stdio.h>
#include <stdlib.h>

#include "agdeblend.h"

int main(void) {
  size_t n_traces = 128;
  size_t n_times[2] = {512, 768};
  int n_patches = 2;
  int volumes[2] = {0, 1};
  int n_dims[2] = {2, 2};
  int **window_shapes = (int **)malloc(2 * sizeof(int *));
  int **coords = (int **)malloc(2 * sizeof(int *));
  int **shapes = (int **)malloc(2 * sizeof(int *));
  int long **shottimes = (int long **)malloc(2 * sizeof(int long *));
  int **channels = (int **)malloc(2 * sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(2 * sizeof(enum AGDTraceType *));
  float **data = (float **)malloc(2 * sizeof(float *));
  float initial_factor = 1.0f;
  int n_its = 2500;
  int print_freq = -1;
  int volume_idx;
  char filename[256];
  FILE *fid;

  for (volume_idx = 0; volume_idx < 2; volume_idx++) {
    window_shapes[volume_idx] = (int *)malloc(n_dims[volume_idx] * sizeof(int));
    coords[volume_idx] = (int *)malloc((n_dims[volume_idx] - 1) * sizeof(int));
    shapes[volume_idx] = (int *)malloc(n_dims[volume_idx] * sizeof(int));
    shottimes[volume_idx] = (int long *)malloc(n_traces * sizeof(int long));
    channels[volume_idx] = (int *)malloc(n_traces * sizeof(int));
    trace_types[volume_idx] =
        (enum AGDTraceType *)malloc(n_traces * sizeof(enum AGDTraceType));
    data[volume_idx] =
        (float *)malloc(n_traces * n_times[volume_idx] * sizeof(float));

    window_shapes[volume_idx][0] = 64;
    window_shapes[volume_idx][1] = 256;
    coords[volume_idx][0] = 0;
    shapes[volume_idx][0] = n_traces;
    shapes[volume_idx][1] = n_times[volume_idx];

    /* Load from input files */
    sprintf(filename, "data_3_shottimes_%d.bin", volume_idx);
    fid = fopen(filename, "rb");
    fread(shottimes[volume_idx], sizeof(long int), n_traces, fid);
    fclose(fid);

    sprintf(filename, "data_3_channels_%d.bin", volume_idx);
    fid = fopen(filename, "rb");
    fread(channels[volume_idx], sizeof(int), n_traces, fid);
    fclose(fid);

    sprintf(filename, "data_3_trace_types_%d.bin", volume_idx);
    fid = fopen(filename, "rb");
    fread(trace_types[volume_idx], sizeof(enum AGDTraceType), n_traces, fid);
    fclose(fid);

    sprintf(filename, "data_3_blended_data_%d.bin", volume_idx);
    fid = fopen(filename, "rb");
    fread(data[volume_idx], sizeof(float), n_traces * n_times[volume_idx], fid);
    fclose(fid);
  }

  /* Deblend */
  agd_deblend(n_patches, volumes, n_dims, (int const *const *)window_shapes,
              (int const *const *)coords, (int const *const *)shapes,
              (long int const *const *)shottimes, (int const *const *)channels,
              (enum AGDTraceType const *const *)trace_types, NULL, NULL, NULL,
              initial_factor, n_its, print_freq, data);

  for (volume_idx = 0; volume_idx < 2; volume_idx++) {
    /* Write to file */
    sprintf(filename, "out/data_3_deblended_data_%d.bin", volume_idx);
    fid = fopen(filename, "wb");
    fwrite(data[volume_idx], sizeof(float), n_traces * n_times[volume_idx],
           fid);
    fclose(fid);

    free(window_shapes[volume_idx]);
    free(coords[volume_idx]);
    free(shapes[volume_idx]);
    free(shottimes[volume_idx]);
    free(channels[volume_idx]);
    free(trace_types[volume_idx]);
    free(data[volume_idx]);
  }
  free(window_shapes);
  free(coords);
  free(shapes);
  free(shottimes);
  free(channels);
  free(trace_types);
  free(data);

  return 0;
}
