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

#include <algorithm>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"
#include "agdeblend.c"
#include "blend_test.cpp"
#include "fk_test.cpp"
#include "model_test.cpp"
#include "wavelet_test.cpp"

TEST(AGDBlend, TwoVolumesOneChannelSameOut) {
  int n_patches = 2;
  int n_patches_out = n_patches;
  int n_traces_per_patch[] = {1, 1};
  int n_traces_per_patch_out[] = {1, 1};
  int n_times_per_patch[] = {3, 3};
  int n_times_per_patch_out[] = {3, 3};
  long int **shottimes = {};
  long int **shottimes_out = {};
  int **channels = {};
  int **channels_out = {};
  enum AGDTraceType **trace_types = {};
  enum AGDTraceType **trace_types_out = {};
  enum AGDBlendMode blend_mode = AGDBlendSum;
  int taper_length = 0;
  AGD_TYPE **data = {};
  AGD_TYPE **data_out = {};
  int n_times = 3;

  shottimes = (long int **)malloc(2 * sizeof(long int *));
  channels = (int **)malloc(2 * sizeof(int *));
  trace_types = (AGDTraceType **)malloc(2 * sizeof(AGDTraceType *));
  data = (AGD_TYPE **)malloc(2 * sizeof(AGD_TYPE *));
  data_out = (AGD_TYPE **)malloc(2 * sizeof(AGD_TYPE *));
  shottimes[0] = (long int *)malloc(sizeof(long int));
  shottimes[1] = (long int *)malloc(sizeof(long int));
  channels[0] = (int *)malloc(sizeof(int));
  channels[1] = (int *)malloc(sizeof(int));
  trace_types[0] = (AGDTraceType *)malloc(sizeof(AGDTraceType));
  trace_types[1] = (AGDTraceType *)malloc(sizeof(AGDTraceType));
  data[0] = (AGD_TYPE *)malloc(n_times * sizeof(AGD_TYPE));
  data[1] = (AGD_TYPE *)malloc(n_times * sizeof(AGD_TYPE));
  data_out[0] = (AGD_TYPE *)malloc(n_times * sizeof(AGD_TYPE));
  data_out[1] = (AGD_TYPE *)malloc(n_times * sizeof(AGD_TYPE));

  shottimes[0][0] = 25;
  shottimes[1][0] = 27;
  channels[0][0] = 123;
  channels[1][0] = 123;
  trace_types[0][0] = AGDLive;
  trace_types[1][0] = AGDLive;
  data[0][0] = 1;
  data[0][1] = 2;
  data[0][2] = 3;
  data[1][0] = 4;
  data[1][1] = 5;
  data[1][2] = 6;

  shottimes_out = shottimes;
  channels_out = channels;
  trace_types_out = trace_types;

  EXPECT_EQ(
      agd_blend(n_patches, n_traces_per_patch, n_times_per_patch, shottimes,
                channels, trace_types, data, blend_mode, taper_length,
                n_patches_out, n_traces_per_patch_out, n_times_per_patch_out,
                shottimes_out, channels_out, trace_types_out, data_out),
      0);
  EXPECT_EQ(data_out[0][0], 1);
  EXPECT_EQ(data_out[0][1], 2);
  EXPECT_EQ(data_out[0][2], 7);
  EXPECT_EQ(data_out[1][0], 7);
  EXPECT_EQ(data_out[1][1], 5);
  EXPECT_EQ(data_out[1][2], 6);

  free(shottimes[0]);
  free(shottimes[1]);
  free(shottimes);
  free(channels[0]);
  free(channels[1]);
  free(channels);
  free(trace_types[0]);
  free(trace_types[1]);
  free(trace_types);
  free(data[0]);
  free(data[1]);
  free(data);
  free(data_out[0]);
  free(data_out[1]);
  free(data_out);
}

TEST(AGDBlend, AllMissing) {
  int n_patches = 1;
  int n_patches_out = n_patches;
  int n_traces_per_patch[] = {1};
  int n_traces_per_patch_out[] = {1};
  int n_times_per_patch[] = {3};
  int n_times_per_patch_out[] = {3};
  long int **shottimes = {};
  long int **shottimes_out = {};
  int **channels = {};
  int **channels_out = {};
  enum AGDTraceType **trace_types = {};
  enum AGDTraceType **trace_types_out = {};
  enum AGDBlendMode blend_mode = AGDBlendSum;
  int taper_length = 0;
  AGD_TYPE **data = {};
  AGD_TYPE **data_out = {};

  shottimes = (long int **)malloc(1 * sizeof(long int *));
  channels = (int **)malloc(1 * sizeof(int *));
  trace_types = (AGDTraceType **)malloc(1 * sizeof(AGDTraceType *));
  data = (AGD_TYPE **)malloc(1 * sizeof(AGD_TYPE *));
  data_out = (AGD_TYPE **)malloc(1 * sizeof(AGD_TYPE *));
  shottimes[0] = (long int *)malloc(sizeof(long int));
  channels[0] = (int *)malloc(sizeof(int));
  trace_types[0] = (AGDTraceType *)malloc(sizeof(AGDTraceType));
  data[0] = (AGD_TYPE *)malloc(3 * sizeof(AGD_TYPE));
  data_out[0] = (AGD_TYPE *)malloc(3 * sizeof(AGD_TYPE));

  shottimes[0][0] = 25;
  channels[0][0] = 123;
  trace_types[0][0] = AGDMissing;
  for (int i = 0; i < 3; i++) {
    data[0][i] = i + 1;
  }

  shottimes_out = shottimes;
  channels_out = channels;
  trace_types_out = trace_types;

  EXPECT_EQ(
      agd_blend(n_patches, n_traces_per_patch, n_times_per_patch, shottimes,
                channels, trace_types, data, blend_mode, taper_length,
                n_patches_out, n_traces_per_patch_out, n_times_per_patch_out,
                shottimes_out, channels_out, trace_types_out, data_out),
      0);
  EXPECT_EQ(data_out[0][0], 0);
  EXPECT_EQ(data_out[0][1], 0);
  EXPECT_EQ(data_out[0][2], 0);

  free(shottimes[0]);
  free(shottimes);
  free(channels[0]);
  free(channels);
  free(trace_types[0]);
  free(trace_types);
  free(data[0]);
  free(data);
  free(data_out[0]);
  free(data_out);
}

TEST(AGDBlend, AllInMissingOutLive) {
  int n_patches = 1;
  int n_patches_out = n_patches;
  int n_traces_per_patch[] = {1};
  int n_traces_per_patch_out[] = {1};
  int n_times_per_patch[] = {3};
  int n_times_per_patch_out[] = {3};
  long int **shottimes = {};
  long int **shottimes_out = {};
  int **channels = {};
  int **channels_out = {};
  enum AGDTraceType **trace_types = {};
  enum AGDTraceType **trace_types_out = {};
  enum AGDBlendMode blend_mode = AGDBlendSum;
  int taper_length = 0;
  AGD_TYPE **data = {};
  AGD_TYPE **data_out = {};

  shottimes = (long int **)malloc(1 * sizeof(long int *));
  channels = (int **)malloc(1 * sizeof(int *));
  trace_types = (AGDTraceType **)malloc(1 * sizeof(AGDTraceType *));
  trace_types_out = (AGDTraceType **)malloc(1 * sizeof(AGDTraceType *));
  data = (AGD_TYPE **)malloc(1 * sizeof(AGD_TYPE *));
  data_out = (AGD_TYPE **)malloc(1 * sizeof(AGD_TYPE *));
  shottimes[0] = (long int *)malloc(sizeof(long int));
  channels[0] = (int *)malloc(sizeof(int));
  trace_types[0] = (AGDTraceType *)malloc(sizeof(AGDTraceType));
  trace_types_out[0] = (AGDTraceType *)malloc(sizeof(AGDTraceType));
  data[0] = (AGD_TYPE *)malloc(3 * sizeof(AGD_TYPE));
  data_out[0] = (AGD_TYPE *)malloc(3 * sizeof(AGD_TYPE));

  shottimes[0][0] = 25;
  channels[0][0] = 123;
  trace_types[0][0] = AGDMissing;
  trace_types_out[0][0] = AGDLive;
  for (int i = 0; i < 3; i++) {
    data[0][i] = i + 1;
  }

  shottimes_out = shottimes;
  channels_out = channels;

  EXPECT_EQ(
      agd_blend(n_patches, n_traces_per_patch, n_times_per_patch, shottimes,
                channels, trace_types, data, blend_mode, taper_length,
                n_patches_out, n_traces_per_patch_out, n_times_per_patch_out,
                shottimes_out, channels_out, trace_types_out, data_out),
      0);
  EXPECT_EQ(data_out[0][0], 0);
  EXPECT_EQ(data_out[0][1], 0);
  EXPECT_EQ(data_out[0][2], 0);

  free(shottimes[0]);
  free(shottimes);
  free(channels[0]);
  free(channels);
  free(trace_types[0]);
  free(trace_types_out[0]);
  free(trace_types);
  free(trace_types_out);
  free(data[0]);
  free(data);
  free(data_out[0]);
  free(data_out);
}

TEST(AGDBlend, AllInLiveOutMissing) {
  int n_patches = 1;
  int n_patches_out = n_patches;
  int n_traces_per_patch[] = {1};
  int n_traces_per_patch_out[] = {1};
  int n_times_per_patch[] = {3};
  int n_times_per_patch_out[] = {3};
  long int **shottimes = {};
  long int **shottimes_out = {};
  int **channels = {};
  int **channels_out = {};
  enum AGDTraceType **trace_types = {};
  enum AGDTraceType **trace_types_out = {};
  enum AGDBlendMode blend_mode = AGDBlendSum;
  int taper_length = 0;
  AGD_TYPE **data = {};
  AGD_TYPE **data_out = {};

  shottimes = (long int **)malloc(1 * sizeof(long int *));
  channels = (int **)malloc(1 * sizeof(int *));
  trace_types = (AGDTraceType **)malloc(1 * sizeof(AGDTraceType *));
  trace_types_out = (AGDTraceType **)malloc(1 * sizeof(AGDTraceType *));
  data = (AGD_TYPE **)malloc(1 * sizeof(AGD_TYPE *));
  data_out = (AGD_TYPE **)malloc(1 * sizeof(AGD_TYPE *));
  shottimes[0] = (long int *)malloc(sizeof(long int));
  channels[0] = (int *)malloc(sizeof(int));
  trace_types[0] = (AGDTraceType *)malloc(sizeof(AGDTraceType));
  trace_types_out[0] = (AGDTraceType *)malloc(sizeof(AGDTraceType));
  data[0] = (AGD_TYPE *)malloc(3 * sizeof(AGD_TYPE));
  data_out[0] = (AGD_TYPE *)malloc(3 * sizeof(AGD_TYPE));

  shottimes[0][0] = 25;
  channels[0][0] = 123;
  trace_types[0][0] = AGDLive;
  trace_types_out[0][0] = AGDMissing;
  for (int i = 0; i < 3; i++) {
    data[0][i] = i + 1;
  }

  shottimes_out = shottimes;
  channels_out = channels;

  EXPECT_EQ(
      agd_blend(n_patches, n_traces_per_patch, n_times_per_patch, shottimes,
                channels, trace_types, data, blend_mode, taper_length,
                n_patches_out, n_traces_per_patch_out, n_times_per_patch_out,
                shottimes_out, channels_out, trace_types_out, data_out),
      0);
  EXPECT_EQ(data_out[0][0], 0);
  EXPECT_EQ(data_out[0][1], 0);
  EXPECT_EQ(data_out[0][2], 0);

  free(shottimes[0]);
  free(shottimes);
  free(channels[0]);
  free(channels);
  free(trace_types[0]);
  free(trace_types_out[0]);
  free(trace_types);
  free(trace_types_out);
  free(data[0]);
  free(data);
  free(data_out[0]);
  free(data_out);
}

TEST(AGDBlend, OneVolumeToTwoWithMute) {
  int n_patches = 1;
  int n_patches_out = 2;
  int n_traces_per_patch[] = {5};
  int n_traces_per_patch_out[] = {3, 2};
  int n_times_per_patch[] = {3};
  int n_times_per_patch_out[] = {2, 3};
  long int **shottimes = {};
  long int **shottimes_out = {};
  int **channels = {};
  int **channels_out = {};
  enum AGDTraceType **trace_types = {};
  enum AGDTraceType **trace_types_out = {};
  enum AGDBlendMode blend_mode = AGDBlendSum;
  int taper_length = 0;
  AGD_TYPE **data = {};
  AGD_TYPE **data_out = {};

  shottimes = (long int **)malloc(1 * sizeof(long int *));
  shottimes_out = (long int **)malloc(2 * sizeof(long int *));
  channels = (int **)malloc(1 * sizeof(int *));
  channels_out = (int **)malloc(2 * sizeof(int *));
  trace_types = (AGDTraceType **)malloc(1 * sizeof(AGDTraceType *));
  trace_types_out = (AGDTraceType **)malloc(2 * sizeof(AGDTraceType *));
  data = (AGD_TYPE **)malloc(1 * sizeof(AGD_TYPE *));
  data_out = (AGD_TYPE **)malloc(2 * sizeof(AGD_TYPE *));
  shottimes[0] = (long int *)malloc(5 * sizeof(long int));
  shottimes_out[0] = (long int *)malloc(3 * sizeof(long int));
  shottimes_out[1] = (long int *)malloc(2 * sizeof(long int));
  channels[0] = (int *)malloc(5 * sizeof(int));
  channels_out[0] = (int *)malloc(3 * sizeof(int));
  channels_out[1] = (int *)malloc(2 * sizeof(int));
  trace_types[0] = (AGDTraceType *)malloc(5 * sizeof(AGDTraceType));
  trace_types_out[0] = (AGDTraceType *)malloc(3 * sizeof(AGDTraceType));
  trace_types_out[1] = (AGDTraceType *)malloc(2 * sizeof(AGDTraceType));
  data[0] = (AGD_TYPE *)malloc(5 * 3 * sizeof(AGD_TYPE));
  data_out[0] = (AGD_TYPE *)malloc(3 * 2 * sizeof(AGD_TYPE));
  data_out[1] = (AGD_TYPE *)malloc(2 * 3 * sizeof(AGD_TYPE));

  shottimes[0][0] = 25;
  shottimes[0][1] = 27;
  shottimes[0][2] = 29;
  shottimes[0][3] = 31;
  shottimes[0][4] = -10;
  shottimes_out[0][0] = 24;
  shottimes_out[0][1] = 26;
  shottimes_out[0][2] = -10;
  shottimes_out[1][0] = 28;
  shottimes_out[1][1] = 30;
  channels[0][0] = 123;
  channels[0][1] = 123;
  channels[0][2] = 123;
  channels[0][3] = 123;
  channels[0][4] = -789;
  channels_out[0][0] = 123;
  channels_out[0][1] = 123;
  channels_out[0][2] = -789;
  channels_out[1][0] = 123;
  channels_out[1][1] = 123;
  trace_types[0][0] = AGDMissing;
  trace_types[0][1] = AGDLive;
  trace_types[0][2] = AGDLive;
  trace_types[0][3] = AGDBad;
  trace_types[0][4] = AGDLive;
  trace_types_out[0][0] = AGDMissing;
  trace_types_out[0][1] = AGDLive;
  trace_types_out[0][2] = AGDLive;
  trace_types_out[1][0] = AGDLive;
  trace_types_out[1][1] = AGDLive;
  for (int i = 0; i < 5 * 3; i++) {
    data[0][i] = i + 1;
  }

  EXPECT_EQ(
      agd_blend(n_patches, n_traces_per_patch, n_times_per_patch, shottimes,
                channels, trace_types, data, blend_mode, taper_length,
                n_patches_out, n_traces_per_patch_out, n_times_per_patch_out,
                shottimes_out, channels_out, trace_types_out, data_out),
      0);

  /* Channel 123:
   * 1  2  3 (Missing, so not added)
   *       4  5  6
   *             7  8  9
   *                   10  11  12 (Bad, so muted)
   * ----------------------------
   * 0  0  4  5  13 8  0   0   0
   *
   * Channel -789:
   * 13 14 15
   */
  EXPECT_EQ(data_out[0][0 * 2 + 0], 0);
  EXPECT_EQ(data_out[0][0 * 2 + 1], 0);
  EXPECT_EQ(data_out[0][1 * 2 + 0], 0);
  EXPECT_EQ(data_out[0][1 * 2 + 1], 4);
  EXPECT_EQ(data_out[0][2 * 2 + 0], 13);
  EXPECT_EQ(data_out[0][2 * 2 + 1], 14);
  EXPECT_EQ(data_out[1][0 * 3 + 0], 5);
  EXPECT_EQ(data_out[1][0 * 3 + 1], 13);
  EXPECT_EQ(data_out[1][0 * 3 + 2], 8);
  EXPECT_EQ(data_out[1][1 * 3 + 0], 8);
  EXPECT_EQ(data_out[1][1 * 3 + 1], 0);
  EXPECT_EQ(data_out[1][1 * 3 + 2], 0);

  free(shottimes[0]);
  free(shottimes_out[0]);
  free(shottimes_out[1]);
  free(shottimes);
  free(shottimes_out);
  free(channels[0]);
  free(channels_out[0]);
  free(channels_out[1]);
  free(channels);
  free(channels_out);
  free(trace_types[0]);
  free(trace_types_out[0]);
  free(trace_types_out[1]);
  free(trace_types);
  free(trace_types_out);
  free(data[0]);
  free(data);
  free(data_out[0]);
  free(data_out[1]);
  free(data_out);
}

TEST(AGDeblend, Deblend) {
  srand(1);
  for (int rand_repeat = 0; rand_repeat < 5; rand_repeat++) {
    int n_patches = 1 + rand() % 3;
    std::vector<int> n_dims(n_patches);
    std::vector<int> volumes(n_patches);
    std::generate(n_dims.begin(), n_dims.end(), [] { return 2 + rand() % 3; });
    std::vector<std::vector<int>> data_shapes(n_patches);
    std::vector<int *> data_shapes_ptr(n_patches, nullptr);
    std::vector<std::vector<int>> window_shapes(n_patches);
    std::vector<int *> window_shapes_ptr(n_patches, nullptr);
    std::vector<std::vector<int>> coords(n_patches);
    std::vector<int *> coords_ptr(n_patches, nullptr);
    std::vector<std::vector<long int>> shottimes(n_patches);
    std::vector<long int *> shottimes_ptr(n_patches, nullptr);
    std::vector<std::vector<int>> channels(n_patches);
    std::vector<int *> channels_ptr(n_patches, nullptr);
    std::vector<std::vector<enum AGDTraceType>> trace_types(n_patches);
    std::vector<enum AGDTraceType *> trace_types_ptr(n_patches, nullptr);
    std::vector<std::vector<AGD_TYPE>> true_data(n_patches);
    std::vector<AGD_TYPE *> true_data_ptr(n_patches, nullptr);
    std::vector<std::vector<AGD_TYPE>> data(n_patches);
    std::vector<AGD_TYPE *> data_ptr(n_patches, nullptr);
    int *wavelet_lengths = nullptr;
    int **wavelet_idxs = nullptr;
    AGD_TYPE **wavelets = nullptr;
    std::vector<int> n_traces(n_patches);
    std::vector<int> n_times(n_patches);
    int n_times_total = 0;

    /* Set data_shapes, window_shapes, n_traces, and n_times */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      data_shapes[patch_idx].resize(n_dims[patch_idx]);
      window_shapes[patch_idx].resize(n_dims[patch_idx]);
      coords[patch_idx].resize(n_dims[patch_idx] - 1);
      for (int dim_idx = 0; dim_idx < n_dims[patch_idx]; dim_idx++) {
        data_shapes[patch_idx][dim_idx] = 8 + rand() % 9;
        window_shapes[patch_idx][dim_idx] =
            8 + rand() % (data_shapes[patch_idx][dim_idx] - 8 + 1);
        /* round down to multiple of 4 if not equal to data shape or time dim */
        if (window_shapes[patch_idx][dim_idx] !=
                data_shapes[patch_idx][dim_idx] ||
            dim_idx == n_dims[patch_idx] - 1) {
          window_shapes[patch_idx][dim_idx] =
              (window_shapes[patch_idx][dim_idx] / 4) * 4;
        }
      }
      n_traces[patch_idx] = 1;
      for (int dim_idx = 0; dim_idx < n_dims[patch_idx] - 1; dim_idx++) {
        n_traces[patch_idx] *= data_shapes[patch_idx][dim_idx];
      }
      n_times[patch_idx] = data_shapes[patch_idx][n_dims[patch_idx] - 1];
      n_times_total += n_traces[patch_idx] * n_times[patch_idx];
    }

    /* Set shottimes, channels, trace_types */
    int shottime_max = n_times_total / 3;
    int channel_max =
        1 + std::accumulate(n_traces.cbegin(), n_traces.cend(), 0) / 10;
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      shottimes[patch_idx].resize(n_traces[patch_idx]);
      channels[patch_idx].resize(n_traces[patch_idx]);
      trace_types[patch_idx].resize(n_traces[patch_idx]);

      std::generate(shottimes[patch_idx].begin(), shottimes[patch_idx].end(),
                    [shottime_max]() { return rand() % shottime_max; });
      std::generate(channels[patch_idx].begin(), channels[patch_idx].end(),
                    [channel_max]() { return rand() % channel_max; });
    }

    /* Set data */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      true_data[patch_idx].resize(n_traces[patch_idx] * n_times[patch_idx]);
      data[patch_idx].resize(n_traces[patch_idx] * n_times[patch_idx]);
      for (int time_idx = 0; time_idx < n_times[patch_idx]; time_idx++) {
        AGD_TYPE time_value = (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX - AGD_HALF;
        for (int trace_idx = 0; trace_idx < n_traces[patch_idx]; trace_idx++) {
          true_data[patch_idx][trace_idx * n_times[patch_idx] + time_idx] =
              time_value;
        }
      }
    }

    std::iota(volumes.begin(), volumes.end(), 0);

    /* Get pointers */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      data_shapes_ptr[patch_idx] = &data_shapes[patch_idx][0];
      window_shapes_ptr[patch_idx] = &window_shapes[patch_idx][0];
      coords_ptr[patch_idx] = &coords[patch_idx][0];
      shottimes_ptr[patch_idx] = &shottimes[patch_idx][0];
      channels_ptr[patch_idx] = &channels[patch_idx][0];
      trace_types_ptr[patch_idx] = &trace_types[patch_idx][0];
      true_data_ptr[patch_idx] = &true_data[patch_idx][0];
      data_ptr[patch_idx] = &data[patch_idx][0];
    }

    /* Blend */
    EXPECT_EQ(
        agd_blend(n_patches, n_traces.data(), n_times.data(),
                  shottimes_ptr.data(), channels_ptr.data(),
                  trace_types_ptr.data(), true_data_ptr.data(), AGDBlendSum, 0,
                  n_patches, n_traces.data(), n_times.data(),
                  shottimes_ptr.data(), channels_ptr.data(),
                  trace_types_ptr.data(), data_ptr.data()),
        0);

    /* Deblend */
    AGD_TYPE initial_factor = AGD_ONE;
    int n_its = 800;
    EXPECT_EQ(agd_deblend(n_patches, volumes.data(), n_dims.data(),
                          window_shapes_ptr.data(), coords_ptr.data(),
                          data_shapes_ptr.data(), shottimes_ptr.data(),
                          channels_ptr.data(), trace_types_ptr.data(),
                          wavelet_lengths, wavelet_idxs, wavelets,
                          initial_factor, n_its, 0, data_ptr.data()),
              0);

    AGD_TYPE error_sq = AGD_ZERO;
    AGD_TYPE true_sq = AGD_ZERO;
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      for (int idx = 0; idx < n_traces[patch_idx] * n_times[patch_idx]; idx++) {
        AGD_TYPE error = data[patch_idx][idx] - true_data[patch_idx][idx];
        error_sq += error * error;
        true_sq += true_data[patch_idx][idx] * true_data[patch_idx][idx];
      }
    }
    EXPECT_EQ(error_sq < true_sq / 10, 1);
  }
}
