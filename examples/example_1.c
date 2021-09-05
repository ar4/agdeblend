#include <stdio.h>
#include <stdlib.h>

#include "agdeblend.h"

int main(void) {
  int const n_patches = 1;
  int const n_traces[1] = {128};
  int const n_times[1] = {512};
  int long **shottimes = (int long **)malloc(sizeof(int long *));
  int **channels = (int **)malloc(sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(sizeof(enum AGDTraceType *));
  float **data = (float **)malloc(sizeof(float *));
  enum AGDBlendMode blend_mode = AGDBlendSum;
  int taper_length = 0;
  FILE *fid;

  shottimes[0] = (int long *)malloc((size_t)n_traces[0] * sizeof(int long));
  channels[0] = (int *)malloc((size_t)n_traces[0] * sizeof(int));
  trace_types[0] = (enum AGDTraceType *)malloc((size_t)n_traces[0] *
                                               sizeof(enum AGDTraceType));
  data[0] =
      (float *)malloc((size_t)n_traces[0] * (size_t)n_times[0] * sizeof(float));

  /* Load from input files */
  fid = fopen("data_1_shottimes.bin", "rb");
  fread(shottimes[0], sizeof(long int), (size_t)n_traces[0], fid);
  fclose(fid);

  fid = fopen("data_1_channels.bin", "rb");
  fread(channels[0], sizeof(int), (size_t)n_traces[0], fid);
  fclose(fid);

  fid = fopen("data_1_trace_types.bin", "rb");
  fread(trace_types[0], sizeof(enum AGDTraceType), (size_t)n_traces[0], fid);
  fclose(fid);

  fid = fopen("data_1_true_data.bin", "rb");
  fread(data[0], sizeof(float), (size_t)n_traces[0] * (size_t)n_times[0], fid);
  fclose(fid);

  /* Blend */
  agd_blend(n_patches, n_traces, n_times, (long int const *const *)shottimes,
            (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types, (float *const *)data,
            blend_mode, taper_length, n_patches, n_traces, n_times,
            (long int const *const *)shottimes, (int const *const *)channels,
            (enum AGDTraceType const *const *)trace_types,
            (float *const *)data);

  /* Write to file */
  fid = fopen("out/data_1_blended_data.bin", "wb");
  fwrite(data[0], sizeof(float), (size_t)n_traces[0] * (size_t)n_times[0], fid);
  fclose(fid);

  free(shottimes[0]);
  free(channels[0]);
  free(trace_types[0]);
  free(data[0]);
  free(shottimes);
  free(channels);
  free(trace_types);
  free(data);

  return 0;
}
