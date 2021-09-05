#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "agdeblend.h"

static void add_overlap(float **data, int n_traces, int n_times, int n_overlap,
                        MPI_Comm comm) {
  int comm_rank;

  MPI_Comm_rank(comm, &comm_rank);

  if (comm_rank == 0) {
    float *buffer = malloc(n_overlap * n_times * sizeof(float));
    int trace_idx;

    MPI_Recv(buffer, n_overlap * n_times, MPI_FLOAT, 1, 0, comm,
             MPI_STATUS_IGNORE);

    for (trace_idx = 0; trace_idx < n_overlap; trace_idx++) {
      int trace0_idx = n_traces - n_overlap + trace_idx;
      int time_idx;
      for (time_idx = 0; time_idx < n_times; time_idx++) {
        data[0][trace0_idx * n_times + time_idx] +=
            buffer[trace_idx * n_times + time_idx];
      }
    }

    free(buffer);
  } else {
    MPI_Send(data[0], n_overlap * n_times, MPI_FLOAT, 0, 0, comm);
  }
}

int main(void) {
  size_t n_traces0 = 80;
  size_t n_traces1 = 64;
  size_t n_times = 512;
  int n_patches = 1;
  int volumes[1] = {0};
  int n_dims[1] = {2};
  int **window_shapes = (int **)malloc(sizeof(int *));
  int **coords = (int **)malloc(sizeof(int *));
  int **shapes = (int **)malloc(sizeof(int *));
  int long **shottimes = (int long **)malloc(sizeof(int long *));
  int **channels = (int **)malloc(sizeof(int *));
  enum AGDTraceType **trace_types =
      (enum AGDTraceType **)malloc(sizeof(enum AGDTraceType *));
  float **data = (float **)malloc(sizeof(float *));
  float initial_factor = 1.0f;
  int n_its = 1000;
  int print_freq = -1;
  int n_overlap;
  int comm_rank;
  int n_traces;
  FILE *fid;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  if (comm_rank == 0) {
    n_traces = n_traces0;
  } else {
    n_traces = n_traces1;
  }

  window_shapes[0] = (int *)malloc(n_dims[0] * sizeof(int));
  coords[0] = (int *)malloc((n_dims[0] - 1) * sizeof(int));
  shapes[0] = (int *)malloc(n_dims[0] * sizeof(int));
  shottimes[0] = (int long *)malloc(n_traces * sizeof(int long));
  channels[0] = (int *)malloc(n_traces * sizeof(int));
  trace_types[0] =
      (enum AGDTraceType *)malloc(n_traces * sizeof(enum AGDTraceType));
  data[0] = (float *)malloc(n_traces * n_times * sizeof(float));

  window_shapes[0][0] = 32;
  window_shapes[0][1] = 256;
  coords[0][0] = comm_rank;
  shapes[0][0] = n_traces;
  shapes[0][1] = n_times;
  n_overlap = window_shapes[0][0] / 2;

  /* Load from input files */
  fid = fopen("data_1_shottimes.bin", "rb");
  fseek(fid, comm_rank * 64 * sizeof(long int), SEEK_SET);
  fread(shottimes[0], sizeof(long int), n_traces, fid);
  fclose(fid);

  fid = fopen("data_1_channels.bin", "rb");
  fseek(fid, comm_rank * 64 * sizeof(int), SEEK_SET);
  fread(channels[0], sizeof(int), n_traces, fid);
  fclose(fid);

  fid = fopen("data_1_trace_types.bin", "rb");
  fseek(fid, comm_rank * 64 * sizeof(int), SEEK_SET);
  fread(trace_types[0], sizeof(enum AGDTraceType), n_traces, fid);
  fclose(fid);

  fid = fopen("data_1_blended_data.bin", "rb");
  fseek(fid, comm_rank * 64 * n_times * sizeof(float), SEEK_SET);
  fread(data[0], sizeof(float), n_traces * n_times, fid);
  fclose(fid);

  /* Deblend */
  agd_deblend(n_patches, volumes, n_dims, (int const *const *)window_shapes,
              (int const *const *)coords, (int const *const *)shapes,
              (long int const *const *)shottimes, (int const *const *)channels,
              (enum AGDTraceType const *const *)trace_types, NULL, NULL, NULL,
              initial_factor, n_its, print_freq, MPI_COMM_WORLD, data);

  /* Write to file */
  /* Add overlap from patch 1 to patch 0 */
  add_overlap(data, n_traces0, n_times, n_overlap, MPI_COMM_WORLD);
  if (comm_rank == 0) {
    fid = fopen("out/data_1_deblended_data.bin", "wb");
    fwrite(data[0], sizeof(float), n_traces0 * n_times, fid);
    fclose(fid);
    MPI_Barrier(MPI_COMM_WORLD);
  } else {
    MPI_Barrier(MPI_COMM_WORLD);
    fid = fopen("out/data_1_deblended_data.bin", "ab");
    fwrite(data[0] + n_overlap * n_times, sizeof(float),
           (n_traces1 - n_overlap) * n_times, fid);
    fclose(fid);
  }

  free(window_shapes[0]);
  free(coords[0]);
  free(shapes[0]);
  free(shottimes[0]);
  free(channels[0]);
  free(trace_types[0]);
  free(data[0]);
  free(window_shapes);
  free(coords);
  free(shapes);
  free(shottimes);
  free(channels);
  free(trace_types);
  free(data);

  MPI_Finalize();

  return 0;
}
