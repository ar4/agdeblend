#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "agdeblend.h"

int main(void) {
  size_t n_traces[7];
  size_t n_times = 512;
  int n_patches = 7;
  int volumes[7] = {0};
  int n_dims[1] = {3};
  int **window_shapes = (int **)malloc(n_patches * sizeof(int *));
  int **coords = (int **)malloc(n_patches * sizeof(int *));
  int **shapes = (int **)malloc(n_patches * sizeof(int *));
  long int **shottimes = (long int **)malloc(n_patches * sizeof(long int *));
  int **channels = (int **)malloc(n_patches * sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(n_patches * sizeof(enum AGDTraceType *));
  float **data = (float **)malloc(n_patches * sizeof(float *));
  float initial_factor = 1.0f;
  int n_its = 1000;
  int print_freq = -1;
  int patch_idx;
  char filename[256];
  FILE *fid;
  int space_shape_all[2] = {32, 40}; /* 40 = 2 * 24 - 8 */
  int n_traces_all = space_shape_all[0] * space_shape_all[1];
  int file_range[2][2][2] = {{{0, 16}, {0, 24}}, {{16, 32}, {16, 40}}};
  int file_shape[2] = {16, 24};
  int n_traces_file = file_shape[0] * file_shape[1];
  int patch_ranges_x[3][2] = {{0, 16}, {8, 24}, {16, 32}};
  int patch_ranges_y[3][2] = {{0, 16}, {8, 32}, {24, 40}};
  int **patch_ranges = malloc(2 * n_patches * sizeof(int *));
  int ix;
  int iy;
  int trace_idx;
  int file_idx;
  long int *shottimes_file =
      (long int *)malloc(n_traces_file * sizeof(long int));
  int *channels_file = (int *)malloc(n_traces_file * sizeof(int));
  enum AGDTraceType *trace_types_file =
      (enum AGDTraceType *)malloc(n_traces_file * sizeof(enum AGDTraceType));
  float *data_file = (float *)malloc(n_traces_file * n_times * sizeof(float));
  long int *shottimes_all = (long int *)malloc(n_traces_all * sizeof(long int));
  int *channels_all = (int *)malloc(n_traces_all * sizeof(int));
  enum AGDTraceType *trace_types_all =
      (enum AGDTraceType *)malloc(n_traces_all * sizeof(enum AGDTraceType));
  float *data_all = (float *)malloc(n_traces_all * n_times * sizeof(float));

  for (trace_idx = 0; trace_idx < n_traces_all; trace_idx++) {
    trace_types_all[trace_idx] = AGDMissing;
  }

  window_shapes[0] = (int *)malloc(n_dims[0] * sizeof(int));
  window_shapes[0][0] = 16;
  window_shapes[0][1] = 16;
  window_shapes[0][2] = 256;

  patch_idx = 0;
  for (ix = 0; ix < 3; ix++) {
    for (iy = 0; iy < 3; iy++) {
      if ((ix == 2 && iy == 0) || (ix == 0 && iy == 2)) continue;
      patch_ranges[2 * patch_idx] = patch_ranges_x[ix];
      patch_ranges[2 * patch_idx + 1] = patch_ranges_y[iy];
      coords[patch_idx] = (int *)malloc((n_dims[0] - 1) * sizeof(int));
      coords[patch_idx][0] = ix;
      coords[patch_idx][1] = iy;
      shapes[patch_idx] = (int *)malloc(n_dims[0] * sizeof(int));
      shapes[patch_idx][0] =
          patch_ranges[2 * patch_idx][1] - patch_ranges[2 * patch_idx][0];
      shapes[patch_idx][1] = patch_ranges[2 * patch_idx + 1][1] -
                             patch_ranges[2 * patch_idx + 1][0];
      shapes[patch_idx][2] = n_times;
      n_traces[patch_idx] = shapes[patch_idx][0] * shapes[patch_idx][1];
      patch_idx++;
    }
  }

  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    shottimes[patch_idx] =
        (long int *)malloc(n_traces[patch_idx] * sizeof(long int));
    channels[patch_idx] = (int *)malloc(n_traces[patch_idx] * sizeof(int));
    trace_types[patch_idx] = (enum AGDTraceType *)malloc(
        n_traces[patch_idx] * sizeof(enum AGDTraceType));
    data[patch_idx] =
        (float *)malloc(n_traces[patch_idx] * n_times * sizeof(float));
  }

  /* Load from input files and add to _all arrays */
  for (file_idx = 0; file_idx < 2; file_idx++) {
    int file_start_x = file_range[file_idx][0][0];
    int file_start_y = file_range[file_idx][1][0];

    sprintf(filename, "data_4_shottimes_%d.bin", file_idx);
    fid = fopen(filename, "rb");
    fread(shottimes_file, sizeof(long int), n_traces_file, fid);
    fclose(fid);

    sprintf(filename, "data_4_channels_%d.bin", file_idx);
    fid = fopen(filename, "rb");
    fread(channels_file, sizeof(int), n_traces_file, fid);
    fclose(fid);

    sprintf(filename, "data_4_trace_types_%d.bin", file_idx);
    fid = fopen(filename, "rb");
    fread(trace_types_file, sizeof(enum AGDTraceType), n_traces_file, fid);
    fclose(fid);

    sprintf(filename, "data_4_blended_data_%d.bin", file_idx);
    fid = fopen(filename, "rb");
    fread(data_file, sizeof(float), n_traces_file * n_times, fid);
    fclose(fid);

    for (ix = 0; ix < file_shape[0]; ix++) {
      int trace_idx_all =
          (file_start_x + ix) * space_shape_all[1] + file_start_y;
      int trace_idx_file = ix * file_shape[1];
      memcpy(shottimes_all + trace_idx_all, shottimes_file + trace_idx_file,
             file_shape[1] * sizeof(long int));
      memcpy(channels_all + trace_idx_all, channels_file + trace_idx_file,
             file_shape[1] * sizeof(int));
      memcpy(trace_types_all + trace_idx_all, trace_types_file + trace_idx_file,
             file_shape[1] * sizeof(enum AGDTraceType));
      memcpy(data_all + trace_idx_all * n_times,
             data_file + trace_idx_file * n_times,
             file_shape[1] * n_times * sizeof(float));
    }
  }

  /* Extract patches from the _all arrays */
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int patch_start_x = patch_ranges[2 * patch_idx][0];
    int patch_start_y = patch_ranges[2 * patch_idx + 1][0];
    for (ix = 0; ix < shapes[patch_idx][0]; ix++) {
      int trace_idx_all =
          (patch_start_x + ix) * space_shape_all[1] + patch_start_y;
      int trace_idx_patch = ix * shapes[patch_idx][1];
      memcpy(shottimes[patch_idx] + trace_idx_patch,
             shottimes_all + trace_idx_all,
             shapes[patch_idx][1] * sizeof(long int));
      memcpy(channels[patch_idx] + trace_idx_patch,
             channels_all + trace_idx_all, shapes[patch_idx][1] * sizeof(int));
      memcpy(trace_types[patch_idx] + trace_idx_patch,
             trace_types_all + trace_idx_all,
             shapes[patch_idx][1] * sizeof(enum AGDTraceType));
      memcpy(data[patch_idx] + trace_idx_patch * n_times,
             data_all + trace_idx_all * n_times,
             shapes[patch_idx][1] * n_times * sizeof(float));
    }
  }

  /* Free the _file and _all arrays */
  free(shottimes_file);
  free(channels_file);
  free(trace_types_file);
  free(data_file);
  free(shottimes_all);
  free(channels_all);
  free(trace_types_all);
  free(data_all);

  /* Deblend */
  agd_deblend(n_patches, volumes, n_dims, (int const *const *)window_shapes,
              (int const *const *)coords, (int const *const *)shapes,
              (long int const *const *)shottimes, (int const *const *)channels,
              (enum AGDTraceType const *const *)trace_types, NULL, NULL, NULL,
              initial_factor, n_its, print_freq, data);

  /* Allocate data_all and add patches to it */
  data_all = (float *)calloc(n_traces_all * n_times, sizeof(float));
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int patch_start_x = patch_ranges[2 * patch_idx][0];
    int patch_start_y = patch_ranges[2 * patch_idx + 1][0];
    for (ix = 0; ix < shapes[patch_idx][0]; ix++) {
      for (iy = 0; iy < shapes[patch_idx][1]; iy++) {
        int trace_idx_all =
            (patch_start_x + ix) * space_shape_all[1] * n_times +
            (patch_start_y + iy) * n_times;
        int trace_idx_patch =
            ix * shapes[patch_idx][1] * n_times + iy * n_times;
        int time_idx;
        for (time_idx = 0; time_idx < n_times; time_idx++) {
          data_all[trace_idx_all + time_idx] +=
              data[patch_idx][trace_idx_patch + time_idx];
        }
      }
    }
  }

  /* Allocate data_file, extract data for each file and write out */
  data_file = (float *)malloc(n_traces_file * n_times * sizeof(float));
  for (file_idx = 0; file_idx < 2; file_idx++) {
    int file_start_x = file_range[file_idx][0][0];
    int file_start_y = file_range[file_idx][1][0];
    for (ix = 0; ix < file_shape[0]; ix++) {
      int trace_idx_all = (file_start_x + ix) * space_shape_all[1] * n_times +
                          file_start_y * n_times;
      int trace_idx_file = ix * file_shape[1] * n_times;
      memcpy(data_file + trace_idx_file, data_all + trace_idx_all,
             file_shape[1] * n_times * sizeof(float));
    }
    sprintf(filename, "out/data_4_deblended_data_%d.bin", file_idx);
    fid = fopen(filename, "wb");
    fwrite(data_file, sizeof(float), n_traces_file * n_times, fid);
    fclose(fid);
  }
  free(data_file);
  free(data_all);

  free(window_shapes[0]);
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    free(coords[patch_idx]);
    free(shapes[patch_idx]);
    free(shottimes[patch_idx]);
    free(channels[patch_idx]);
    free(trace_types[patch_idx]);
    free(data[patch_idx]);
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
