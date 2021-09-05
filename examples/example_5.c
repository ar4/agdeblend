#include <stdio.h>
#include <stdlib.h>

#include "agdeblend.h"

static void add_overlap(float **data, int n_traces, int n_times,
                        int n_overlap) {
  int trace_idx;
  for (trace_idx = 0; trace_idx < n_overlap; trace_idx++) {
    int trace0_idx = n_traces - n_overlap + trace_idx;
    int time_idx;
    for (time_idx = 0; time_idx < n_times; time_idx++) {
      data[0][trace0_idx * n_times + time_idx] +=
          data[1][trace_idx * n_times + time_idx];
    }
  }
}

int main(void) {
  size_t n_traces0 = 80;
  size_t n_traces1 = 64;
  size_t n_times = 512;
  int n_patches = 2;
  int volumes[2] = {0, 0};
  int n_dims[1] = {2};
  int **window_shapes = (int **)malloc(sizeof(int *));
  int **coords = (int **)malloc(2 * sizeof(int *));
  int **shapes = (int **)malloc(2 * sizeof(int *));
  int long **shottimes = (int long **)malloc(2 * sizeof(int long *));
  int **channels = (int **)malloc(2 * sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(2 * sizeof(enum AGDTraceType *));
  float **data = (float **)malloc(2 * sizeof(float *));
  float initial_factor = 1.0f;
  int n_its = 1000;
  int print_freq = -1;
  int n_overlap;
  FILE *fid;

  window_shapes[0] = (int *)malloc(n_dims[0] * sizeof(int));
  coords[0] = (int *)malloc((n_dims[0] - 1) * sizeof(int));
  coords[1] = (int *)malloc((n_dims[0] - 1) * sizeof(int));
  shapes[0] = (int *)malloc(n_dims[0] * sizeof(int));
  shapes[1] = (int *)malloc(n_dims[0] * sizeof(int));
  shottimes[0] = (int long *)malloc(n_traces0 * sizeof(int long));
  shottimes[1] = (int long *)malloc(n_traces1 * sizeof(int long));
  channels[0] = (int *)malloc(n_traces0 * sizeof(int));
  channels[1] = (int *)malloc(n_traces1 * sizeof(int));
  trace_types[0] =
      (enum AGDTraceType *)malloc(n_traces0 * sizeof(enum AGDTraceType));
  trace_types[1] =
      (enum AGDTraceType *)malloc(n_traces1 * sizeof(enum AGDTraceType));
  data[0] = (float *)malloc(n_traces0 * n_times * sizeof(float));
  data[1] = (float *)malloc(n_traces1 * n_times * sizeof(float));

  window_shapes[0][0] = 32;
  window_shapes[0][1] = 256;
  coords[0][0] = 0;
  coords[1][0] = 1;
  shapes[0][0] = n_traces0;
  shapes[1][0] = n_traces1;
  shapes[0][1] = n_times;
  shapes[1][1] = n_times;
  n_overlap = window_shapes[0][0] / 2;

  /* Load from input files */
  fid = fopen("data_1_shottimes.bin", "rb");
  fread(shottimes[0], sizeof(long int), n_traces0, fid);
  fseek(fid, 64 * sizeof(long int), SEEK_SET);
  fread(shottimes[1], sizeof(long int), n_traces1, fid);
  fclose(fid);

  fid = fopen("data_1_channels.bin", "rb");
  fread(channels[0], sizeof(int), n_traces0, fid);
  fseek(fid, 64 * sizeof(int), SEEK_SET);
  fread(channels[1], sizeof(int), n_traces1, fid);
  fclose(fid);

  fid = fopen("data_1_trace_types.bin", "rb");
  fread(trace_types[0], sizeof(enum AGDTraceType), n_traces0, fid);
  fseek(fid, 64 * sizeof(int), SEEK_SET);
  fread(trace_types[1], sizeof(enum AGDTraceType), n_traces1, fid);
  fclose(fid);

  fid = fopen("data_1_blended_data.bin", "rb");
  fread(data[0], sizeof(float), n_traces0 * n_times, fid);
  fseek(fid, 64 * n_times * sizeof(float), SEEK_SET);
  fread(data[1], sizeof(float), n_traces1 * n_times, fid);
  fclose(fid);

  /* Deblend */
  agd_deblend(n_patches, volumes, n_dims, (int const *const *)window_shapes,
              (int const *const *)coords, (int const *const *)shapes,
              (long int const *const *)shottimes, (int const *const *)channels,
              (enum AGDTraceType const *const *)trace_types, NULL, NULL, NULL,
              initial_factor, n_its, print_freq, data);

  /* Write to file */
  /* Add overlap from patch 1 to patch 0 */
  add_overlap(data, n_traces0, n_times, n_overlap);
  fid = fopen("out/data_1_deblended_data.bin", "wb");
  fwrite(data[0], sizeof(float), n_traces0 * n_times, fid);
  fwrite(data[1] + n_overlap * n_times, sizeof(float),
         (n_traces1 - n_overlap) * n_times, fid);
  fclose(fid);

  free(window_shapes[0]);
  free(coords[0]);
  free(coords[1]);
  free(shapes[0]);
  free(shapes[1]);
  free(shottimes[0]);
  free(shottimes[1]);
  free(channels[0]);
  free(channels[1]);
  free(trace_types[0]);
  free(trace_types[1]);
  free(data[0]);
  free(data[1]);
  free(window_shapes);
  free(coords);
  free(shapes);
  free(shottimes);
  free(channels);
  free(trace_types);
  free(data);

  return 0;
}
