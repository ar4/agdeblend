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

#include "agdeblend.c"
#include "gtest/gtest.h"

TEST(Model, SetSpaceWindowsMPI) {
  MPI_Comm newcomm;
  int comm_rank;
  int comm_size;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if (comm_size < 2) {
    fprintf(
        stderr,
        "WARNING: More than 2 ranks are required for this test. SKIPPING!\n");
    return;
  }
  MPI_Comm_split(MPI_COMM_WORLD, comm_rank < 2 ? 0 : MPI_UNDEFINED, 0,
                 &newcomm);
  MPI_Comm_rank(newcomm, &comm_rank);
  MPI_Comm_size(newcomm, &comm_size);

  /* Volume 0:
   * 0 0/1 1/2 2   3 3/4 4/5 5/6 6/7 7
   *
   * Volume 1:
   * 11 11   8 8/9 9/10 10
   * 11 11   8 8/9 9/10 10
   *
   *
   *
   * 12 12
   * 12 12
   */

  if (comm_rank == 0) {
    int const n_patches = 2;
    int const volumes[2] = {0, 1};
    int const n_dims[2] = {2, 3};
    std::vector<std::vector<int>> coords = {{0}, {1, 1}};
    std::vector<int *> coords_ptr = {&coords[0][0], &coords[1][0]};
    std::vector<std::vector<int>> data_shapes = {{4, 8}, {2, 4, 7}};
    std::vector<int *> data_shapes_ptr = {&data_shapes[0][0],
                                          &data_shapes[1][0]};
    std::vector<std::vector<int>> window_shapes = {{2, 8}, {2, 2, 7}};
    std::vector<int *> window_shapes_ptr = {&window_shapes[0][0],
                                            &window_shapes[1][0]};
    int n_windows;
    struct ModelWindow *windows;

    EXPECT_EQ(
        set_space_windows(n_patches, volumes, n_dims, window_shapes_ptr.data(),
                          coords_ptr.data(), data_shapes_ptr.data(), newcomm,
                          &n_windows, &windows),
        0);

    EXPECT_EQ(n_windows, 6);

    struct SpaceWindow window;

    window = windows[0].space_window;
    EXPECT_EQ(window.patch_idx, 0);
    EXPECT_EQ(window.n_dims, 2);
    EXPECT_EQ(window.n_traces, 2);
    EXPECT_EQ(window.coords[0], 0);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.use_taper[0], 0);
    EXPECT_EQ(window.use_taper[1], 1);

    window = windows[1].space_window;
    EXPECT_EQ(window.patch_idx, 0);
    EXPECT_EQ(window.n_dims, 2);
    EXPECT_EQ(window.n_traces, 2);
    EXPECT_EQ(window.coords[0], 1);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.use_taper[0], 1);
    EXPECT_EQ(window.use_taper[1], 1);

    window = windows[2].space_window;
    EXPECT_EQ(window.patch_idx, 0);
    EXPECT_EQ(window.n_dims, 2);
    EXPECT_EQ(window.n_traces, 2);
    EXPECT_EQ(window.coords[0], 2);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.use_taper[0], 1);
    EXPECT_EQ(window.use_taper[1], 1);

    window = windows[3].space_window;
    EXPECT_EQ(window.patch_idx, 1);
    EXPECT_EQ(window.n_dims, 3);
    EXPECT_EQ(window.n_traces, 4);
    EXPECT_EQ(window.coords[0], 0);
    EXPECT_EQ(window.coords[1], 0);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.space_shape[1], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.taper_length[1], 1);
    EXPECT_EQ(window.use_taper[0], 0);
    EXPECT_EQ(window.use_taper[1], 0);
    EXPECT_EQ(window.use_taper[2], 1);
    EXPECT_EQ(window.use_taper[3], 1);

    window = windows[4].space_window;
    EXPECT_EQ(window.patch_idx, 1);
    EXPECT_EQ(window.n_dims, 3);
    EXPECT_EQ(window.n_traces, 4);
    EXPECT_EQ(window.coords[0], 0);
    EXPECT_EQ(window.coords[1], 1);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.space_shape[1], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.taper_length[1], 1);
    EXPECT_EQ(window.use_taper[0], 0);
    EXPECT_EQ(window.use_taper[1], 0);
    EXPECT_EQ(window.use_taper[2], 1);
    EXPECT_EQ(window.use_taper[3], 1);

    window = windows[5].space_window;
    EXPECT_EQ(window.patch_idx, 1);
    EXPECT_EQ(window.n_dims, 3);
    EXPECT_EQ(window.n_traces, 4);
    EXPECT_EQ(window.coords[0], 0);
    EXPECT_EQ(window.coords[1], 2);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.space_shape[1], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.taper_length[1], 1);
    EXPECT_EQ(window.use_taper[0], 0);
    EXPECT_EQ(window.use_taper[1], 0);
    EXPECT_EQ(window.use_taper[2], 1);
    EXPECT_EQ(window.use_taper[3], 0);

    for (int window_idx = 0; window_idx < n_windows; window_idx++) {
      struct ModelWindow *const window = windows + window_idx;
      free_space_window(&(window->space_window));
    }
    free(windows);
  } else if (comm_rank == 1) {
    int const n_patches = 3;
    int const volumes[5] = {0, 1, 1};
    int const n_dims[2] = {2, 3};
    std::vector<std::vector<int>> coords = {{1}, {1, 0}, {3, 1}};
    std::vector<int *> coords_ptr = {&coords[0][0], &coords[1][0],
                                     &coords[2][0]};
    std::vector<std::vector<int>> data_shapes = {{6, 8}, {2, 2, 7}, {2, 2, 7}};
    std::vector<int *> data_shapes_ptr = {
        &data_shapes[0][0], &data_shapes[1][0], &data_shapes[2][0]};
    std::vector<std::vector<int>> window_shapes = {{2, 8}, {2, 2, 7}};
    std::vector<int *> window_shapes_ptr = {&window_shapes[0][0],
                                            &window_shapes[1][0]};
    int n_windows;
    struct ModelWindow *windows;

    EXPECT_EQ(
        set_space_windows(n_patches, volumes, n_dims, window_shapes_ptr.data(),
                          coords_ptr.data(), data_shapes_ptr.data(), newcomm,
                          &n_windows, &windows),
        0);

    EXPECT_EQ(n_windows, 7);

    struct SpaceWindow window;

    window = windows[0].space_window;
    EXPECT_EQ(window.patch_idx, 0);
    EXPECT_EQ(window.n_dims, 2);
    EXPECT_EQ(window.n_traces, 2);
    EXPECT_EQ(window.coords[0], 0);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.use_taper[0], 1);
    EXPECT_EQ(window.use_taper[1], 1);

    window = windows[1].space_window;
    EXPECT_EQ(window.patch_idx, 0);
    EXPECT_EQ(window.n_dims, 2);
    EXPECT_EQ(window.n_traces, 2);
    EXPECT_EQ(window.coords[0], 1);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.use_taper[0], 1);
    EXPECT_EQ(window.use_taper[1], 1);

    window = windows[2].space_window;
    EXPECT_EQ(window.patch_idx, 0);
    EXPECT_EQ(window.n_dims, 2);
    EXPECT_EQ(window.n_traces, 2);
    EXPECT_EQ(window.coords[0], 2);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.use_taper[0], 1);
    EXPECT_EQ(window.use_taper[1], 1);

    window = windows[3].space_window;
    EXPECT_EQ(window.patch_idx, 0);
    EXPECT_EQ(window.n_dims, 2);
    EXPECT_EQ(window.n_traces, 2);
    EXPECT_EQ(window.coords[0], 3);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.use_taper[0], 1);
    EXPECT_EQ(window.use_taper[1], 1);

    window = windows[4].space_window;
    EXPECT_EQ(window.patch_idx, 0);
    EXPECT_EQ(window.n_dims, 2);
    EXPECT_EQ(window.n_traces, 2);
    EXPECT_EQ(window.coords[0], 4);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.use_taper[0], 1);
    EXPECT_EQ(window.use_taper[1], 0);

    window = windows[5].space_window;
    EXPECT_EQ(window.patch_idx, 1);
    EXPECT_EQ(window.n_dims, 3);
    EXPECT_EQ(window.n_traces, 4);
    EXPECT_EQ(window.coords[0], 0);
    EXPECT_EQ(window.coords[1], 0);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.space_shape[1], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.taper_length[1], 1);
    EXPECT_EQ(window.use_taper[0], 0);
    EXPECT_EQ(window.use_taper[1], 0);
    EXPECT_EQ(window.use_taper[2], 0);
    EXPECT_EQ(window.use_taper[3], 1);

    window = windows[6].space_window;
    EXPECT_EQ(window.patch_idx, 2);
    EXPECT_EQ(window.n_dims, 3);
    EXPECT_EQ(window.n_traces, 4);
    EXPECT_EQ(window.coords[0], 0);
    EXPECT_EQ(window.coords[1], 0);
    EXPECT_EQ(window.space_shape[0], 2);
    EXPECT_EQ(window.space_shape[1], 2);
    EXPECT_EQ(window.taper_length[0], 1);
    EXPECT_EQ(window.taper_length[1], 1);
    EXPECT_EQ(window.use_taper[0], 0);
    EXPECT_EQ(window.use_taper[1], 0);
    EXPECT_EQ(window.use_taper[2], 0);
    EXPECT_EQ(window.use_taper[3], 0);

    for (int window_idx = 0; window_idx < n_windows; window_idx++) {
      struct ModelWindow *const window = windows + window_idx;
      free_space_window(&(window->space_window));
    }
    free(windows);
  }

  MPI_Comm_free(&newcomm);
  // MPI_Finalize();
}

TEST(Blend, TwoProcessesMPI) {
  MPI_Comm newcomm;
  int comm_rank;
  int comm_size;

  // MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if (comm_size < 2) {
    fprintf(
        stderr,
        "WARNING: More than 2 ranks are required for this test. SKIPPING!\n");
    return;
  }
  MPI_Comm_split(MPI_COMM_WORLD, comm_rank < 2 ? 0 : MPI_UNDEFINED, 0,
                 &newcomm);
  MPI_Comm_rank(newcomm, &comm_rank);
  MPI_Comm_size(newcomm, &comm_size);

  if (comm_rank == 0) {
    int n_patches = 2;
    int n_traces_per_patch[] = {2, 3};
    int n_times_per_patch[] = {100, 50};
    std::vector<std::vector<long int>> shottimes = {{25, -100}, {35, 45, 55}};
    std::vector<long int *> shottimes_ptr = {&shottimes[0][0],
                                             &shottimes[1][0]};
    std::vector<std::vector<int>> channels = {{5, 71}, {4, 5, 5}};
    std::vector<int *> channels_ptr = {&channels[0][0], &channels[1][0]};
    std::vector<std::vector<enum AGDTraceType>> trace_types = {
        {AGDLive, AGDMissing}, {AGDLive, AGDLive, AGDBad}};
    std::vector<enum AGDTraceType *> trace_types_ptr = {&trace_types[0][0],
                                                        &trace_types[1][0]};
    struct BlendConfig blend_config = {};

    /*
       Channel 4:
       35 - 85 live

       Channel 5:
       25 - 125 (25-55 live, 55-105 bad, 105-125 live)
     */

    EXPECT_EQ(set_blend_config(n_patches, n_traces_per_patch, n_times_per_patch,
                               shottimes_ptr.data(), channels_ptr.data(),
                               trace_types_ptr.data(), newcomm, &blend_config),
              0);
    EXPECT_EQ(blend_config.n_channels, 2);
    EXPECT_EQ(blend_config.unique_channel_values[0], 5);
    EXPECT_EQ(blend_config.unique_channel_values[1], 4);
    EXPECT_EQ(blend_config.channels_intervals[0].n_intervals, 1);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].start, 25);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].stop, 25 + 100);
    EXPECT_EQ(blend_config.channels_intervals[1].n_intervals, 1);
    EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].start, 35);
    EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].stop, 35 + 50);
    EXPECT_EQ(blend_config.n_overlaps[0], 0);
    EXPECT_EQ(blend_config.n_overlaps[1], 2);
    EXPECT_EQ(blend_config.ranks_coords[1][0].channel_idx, 1);
    EXPECT_EQ(blend_config.ranks_coords[1][0].interval_idx, 0);
    EXPECT_EQ(blend_config.ranks_coords[1][0].interval_start, 5);
    EXPECT_EQ(blend_config.ranks_coords[1][0].n_times, 30);
    EXPECT_EQ(blend_config.ranks_coords[1][1].channel_idx, 0);
    EXPECT_EQ(blend_config.ranks_coords[1][1].interval_idx, 0);
    EXPECT_EQ(blend_config.ranks_coords[1][1].interval_start, 0);
    EXPECT_EQ(blend_config.ranks_coords[1][1].n_times, 5);
    EXPECT_EQ(blend_config.comm_size, 2);
    EXPECT_EQ(blend_config.n_mute_coords, 2);
    EXPECT_EQ(blend_config.mute_coords[0].channel_idx, 0);
    EXPECT_EQ(blend_config.mute_coords[0].interval_idx, 0);
    EXPECT_EQ(blend_config.mute_coords[0].interval_start, 30);
    EXPECT_EQ(blend_config.mute_coords[0].n_times, 50);
    EXPECT_EQ(blend_config.mute_coords[1].channel_idx, 1);
    EXPECT_EQ(blend_config.mute_coords[1].interval_idx, 0);
    EXPECT_EQ(blend_config.mute_coords[1].interval_start, 5);
    EXPECT_EQ(blend_config.mute_coords[1].n_times, 10);

    std::vector<std::vector<AGD_TYPE>> data(n_patches);
    data[0].resize(n_traces_per_patch[0] * n_times_per_patch[0]);
    data[1].resize(n_traces_per_patch[1] * n_times_per_patch[1]);
    std::iota(data[0].begin(), data[0].end(), 1);
    std::iota(data[1].begin(), data[1].end(), 500);
    std::vector<AGD_TYPE *> data_ptr = {&data[0][0], &data[1][0]};

    std::vector<struct BlendParams> blend_params(n_patches);
    AGD_TYPE ***blended = nullptr;

    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      EXPECT_EQ(set_blend_params(
                    n_traces_per_patch[patch_idx], n_times_per_patch[patch_idx],
                    shottimes[patch_idx].data(), channels[patch_idx].data(),
                    trace_types[patch_idx].data(), &blend_config,
                    &blend_params[patch_idx]),
                0);
    }

    /* allocate blended */
    EXPECT_EQ(allocate_blended(&blend_config, &blended), 0);
    zero_blended(&blend_config, (AGD_TYPE *const *const *)blended);

    /* blend overwrite */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      blend_overwrite_forward(data[patch_idx].data(),
                              blend_params.data() + patch_idx,
                              (AGD_TYPE *const *const *)blended);
    }
    EXPECT_EQ(blend_overwrite_forward_mpi(&blend_config,
                                          (AGD_TYPE *const *const *)blended),
              0);

    for (int i = 25 - 25; i < 30 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1015 + i));
    }
    for (int i = 30 - 25; i < 45 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1 + i));
    }
    for (int i = 45 - 25; i < 55 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(550 - (45 - 25) + i));
    }
    for (int i = 55 - 25; i < 105 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(600 - (55 - 25) + i));
    }
    for (int i = 105 - 25; i < 100; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1 + i));
    }

    for (int i = 35 - 35; i < 40 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(500 + i));
    }
    for (int i = 40 - 35; i < 50 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(1000 - (40 - 35) + i));
    }
    for (int i = 50 - 35; i < 70 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(1520 - (50 - 35) + i));
    }
    for (int i = 70 - 35; i < 85 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(500 + i));
    }

    /* Overwrite with mute */
    apply_mute(&blend_config, (AGD_TYPE *const *const *)blended);

    for (int i = 25 - 25; i < 30 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1015 + i));
    }
    for (int i = 30 - 25; i < 45 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1 + i));
    }
    for (int i = 45 - 25; i < 55 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(550 - (45 - 25) + i));
    }
    for (int i = 55 - 25; i < 105 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], AGD_ZERO);
    }
    for (int i = 105 - 25; i < 100; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1 + i));
    }

    for (int i = 35 - 35; i < 40 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(500 + i));
    }
    for (int i = 40 - 35; i < 50 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], AGD_ZERO);
    }
    for (int i = 50 - 35; i < 70 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(1520 - (50 - 35) + i));
    }
    for (int i = 70 - 35; i < 85 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(500 + i));
    }

    zero_blended(&blend_config, (AGD_TYPE *const *const *)blended);

    /* blend sum */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      blend_sum_forward(data[patch_idx].data(), blend_params.data() + patch_idx,
                        (AGD_TYPE *const *const *)blended);
    }
    EXPECT_EQ(
        blend_sum_forward_mpi(&blend_config, (AGD_TYPE *const *const *)blended),
        0);

    for (int i = 25 - 25; i < 30 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1015 + 1 + 2 * i));
    }
    for (int i = 30 - 25; i < 45 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1 + i));
    }
    for (int i = 45 - 25; i < 55 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i],
                      (AGD_TYPE)(1 + i + 550 - (45 - 25) + i));
    }
    for (int i = 55 - 25; i < 95 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1 + i + 550 - (45 - 25) + i +
                                                   600 - (55 - 25) + i));
    }
    for (int i = 95 - 25; i < 105 - 25; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i],
                      (AGD_TYPE)(1 + i + 600 - (55 - 25) + i));
    }
    for (int i = 105 - 25; i < 100; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1 + i));
    }

    for (int i = 35 - 35; i < 40 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(500 + i));
    }
    for (int i = 40 - 35; i < 50 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i],
                      (AGD_TYPE)(500 + i + 1000 - (40 - 35) + i));
    }
    for (int i = 50 - 35; i < 70 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i],
                      (AGD_TYPE)(500 + i + 1520 - (50 - 35) + i));
    }
    for (int i = 70 - 35; i < 85 - 35; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(500 + i));
    }

    free_blended(&blend_config, &blended);
    for (auto &blend_param : blend_params) {
      free_blend_params(&blend_param);
    }

    free_blend_config(&blend_config);

  } else if (comm_rank == 1) {
    int n_patches = 2;
    int n_traces_per_patch[] = {2, 2};
    int n_times_per_patch[] = {10, 20};
    std::vector<std::vector<long int>> shottimes = {{40, 20}, {0, 50}};
    std::vector<long int *> shottimes_ptr = {&shottimes[0][0],
                                             &shottimes[1][0]};
    std::vector<std::vector<int>> channels = {{4, 5}, {4, 4}};
    std::vector<int *> channels_ptr = {&channels[0][0], &channels[1][0]};
    std::vector<std::vector<enum AGDTraceType>> trace_types = {
        {AGDBad, AGDLive}, {AGDLive, AGDLive}};
    std::vector<enum AGDTraceType *> trace_types_ptr = {&trace_types[0][0],
                                                        &trace_types[1][0]};
    struct BlendConfig blend_config = {};

    /*
       Channel 4:
       0 - 20 live
       40 - 50 bad
       50 - 70 live

       Channel 5:
       20 - 30 live
     */

    EXPECT_EQ(set_blend_config(n_patches, n_traces_per_patch, n_times_per_patch,
                               shottimes_ptr.data(), channels_ptr.data(),
                               trace_types_ptr.data(), newcomm, &blend_config),
              0);
    EXPECT_EQ(blend_config.n_channels, 2);
    EXPECT_EQ(blend_config.unique_channel_values[0], 4);
    EXPECT_EQ(blend_config.unique_channel_values[1], 5);
    EXPECT_EQ(blend_config.channels_intervals[0].n_intervals, 2);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].start, 40);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].stop, 70);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[1].start, 0);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[1].stop, 20);
    EXPECT_EQ(blend_config.channels_intervals[1].n_intervals, 1);
    EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].start, 20);
    EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].stop, 30);
    EXPECT_EQ(blend_config.n_overlaps[0], 2);
    EXPECT_EQ(blend_config.n_overlaps[1], 0);
    EXPECT_EQ(blend_config.ranks_coords[0][0].channel_idx, 0);
    EXPECT_EQ(blend_config.ranks_coords[0][0].interval_idx, 0);
    EXPECT_EQ(blend_config.ranks_coords[0][0].interval_start, 0);
    EXPECT_EQ(blend_config.ranks_coords[0][0].n_times, 30);
    EXPECT_EQ(blend_config.ranks_coords[0][1].channel_idx, 1);
    EXPECT_EQ(blend_config.ranks_coords[0][1].interval_idx, 0);
    EXPECT_EQ(blend_config.ranks_coords[0][1].interval_start, 5);
    EXPECT_EQ(blend_config.ranks_coords[0][1].n_times, 5);
    EXPECT_EQ(blend_config.comm_size, 2);
    EXPECT_EQ(blend_config.n_mute_coords, 1);
    EXPECT_EQ(blend_config.mute_coords[0].channel_idx, 0);
    EXPECT_EQ(blend_config.mute_coords[0].interval_idx, 0);
    EXPECT_EQ(blend_config.mute_coords[0].interval_start, 0);
    EXPECT_EQ(blend_config.mute_coords[0].n_times, 10);

    std::vector<std::vector<AGD_TYPE>> data(n_patches);
    data[0].resize(n_traces_per_patch[0] * n_times_per_patch[0]);
    data[1].resize(n_traces_per_patch[1] * n_times_per_patch[1]);
    std::iota(data[0].begin(), data[0].end(), 1000);
    std::iota(data[1].begin(), data[1].end(), 1500);
    std::vector<AGD_TYPE *> data_ptr = {&data[0][0], &data[1][0]};

    std::vector<struct BlendParams> blend_params(n_patches);
    AGD_TYPE ***blended = nullptr;

    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      EXPECT_EQ(set_blend_params(
                    n_traces_per_patch[patch_idx], n_times_per_patch[patch_idx],
                    shottimes[patch_idx].data(), channels[patch_idx].data(),
                    trace_types[patch_idx].data(), &blend_config,
                    &blend_params[patch_idx]),
                0);
    }

    /* allocate blended */
    EXPECT_EQ(allocate_blended(&blend_config, &blended), 0);
    zero_blended(&blend_config, (AGD_TYPE *const *const *)blended);

    /* blend overwrite */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      blend_overwrite_forward(data[patch_idx].data(),
                              blend_params.data() + patch_idx,
                              (AGD_TYPE *const *const *)blended);
    }
    EXPECT_EQ(blend_overwrite_forward_mpi(&blend_config,
                                          (AGD_TYPE *const *const *)blended),
              0);

    for (int i = 40 - 40; i < 50 - 40; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1000 + i));
    }
    for (int i = 50 - 40; i < 70 - 40; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1520 - 10 + i));
    }
    for (int i = 0; i < 20; i++) {
      EXPECT_FLOAT_EQ(blended[0][1][i], (AGD_TYPE)(1500 + i));
    }

    for (int i = 20 - 20; i < 30 - 20; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(1010 + i));
    }

    /* Overwrite with mute */
    apply_mute(&blend_config, (AGD_TYPE *const *const *)blended);

    for (int i = 40 - 40; i < 50 - 40; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], AGD_ZERO);
    }
    for (int i = 50 - 40; i < 70 - 40; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(1520 - 10 + i));
    }
    for (int i = 0; i < 20; i++) {
      EXPECT_FLOAT_EQ(blended[0][1][i], (AGD_TYPE)(1500 + i));
    }

    for (int i = 20 - 20; i < 30 - 20; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(1010 + i));
    }

    zero_blended(&blend_config, (AGD_TYPE *const *const *)blended);

    /* blend sum */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      blend_sum_forward(data[patch_idx].data(), blend_params.data() + patch_idx,
                        (AGD_TYPE *const *const *)blended);
    }
    EXPECT_EQ(
        blend_sum_forward_mpi(&blend_config, (AGD_TYPE *const *const *)blended),
        0);

    for (int i = 40 - 40; i < 50 - 40; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(505 + i + 1000 + i));
    }
    for (int i = 50 - 40; i < 70 - 40; i++) {
      EXPECT_FLOAT_EQ(blended[0][0][i], (AGD_TYPE)(505 + i + 1520 - 10 + i));
    }
    for (int i = 0; i < 20; i++) {
      EXPECT_FLOAT_EQ(blended[0][1][i], (AGD_TYPE)(1500 + i));
    }

    for (int i = 20 - 20; i < 25 - 20; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(1010 + i));
    }
    for (int i = 25 - 20; i < 30 - 20; i++) {
      EXPECT_FLOAT_EQ(blended[1][0][i], (AGD_TYPE)(1 + (i - 5) + 1010 + i));
    }

    free_blended(&blend_config, &blended);
    for (auto &blend_param : blend_params) {
      free_blend_params(&blend_param);
    }

    free_blend_config(&blend_config);
  }

  MPI_Comm_free(&newcomm);
  // MPI_Finalize();
}

TEST(Blend, MeanForward) {
  MPI_Comm newcomm;
  int comm_rank;
  int comm_size;

  // MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if (comm_size < 2) {
    fprintf(
        stderr,
        "WARNING: More than 2 ranks are required for this test. SKIPPING!\n");
    return;
  }
  MPI_Comm_split(MPI_COMM_WORLD, comm_rank < 2 ? 0 : MPI_UNDEFINED, 0,
                 &newcomm);
  MPI_Comm_rank(newcomm, &comm_rank);
  MPI_Comm_size(newcomm, &comm_size);

  /* Patch 0: (2 traces, 10 times)
   * Channel 0:
   *  1W0  2W1  3W2  4    5    6    7    8W2  9W1 10W0
   *           11W0 12W1 13W2 14   15   16   17   18W2 19W1 20W0
   *
   * Patch 1: (2 traces, 11 times)
   * Channel 1:
   * 50   51   52   53   54   55   56   57   58   59   60
   * Channel 0:
   *      61W0 62W1 63W2 64   65   66   67   68   69W2 70W1 71W0
   */

  if (comm_rank == 0) {
    int n_patches = 1;
    int n_traces_per_patch[1] = {2};
    int n_times_per_patch[1] = {10};
    int taper_length = 3;
    AGD_TYPE weight[3] = {AGD_HALF * AGD_SQRT(AGD_TWO - AGD_SQRT(AGD_TWO)),
                          AGD_HALF * AGD_SQRT(AGD_TWO),
                          AGD_HALF * AGD_SQRT(AGD_TWO + AGD_SQRT(AGD_TWO))};
    std::vector<std::vector<long int>> shottimes(1);
    std::vector<long int *> shottimes_ptr(n_patches, nullptr);
    std::vector<std::vector<int>> channels(1);
    std::vector<int *> channels_ptr(n_patches, nullptr);
    std::vector<std::vector<enum AGDTraceType>> trace_types(1);
    std::vector<enum AGDTraceType *> trace_types_ptr(n_patches, nullptr);
    std::vector<std::vector<AGD_TYPE>> data(1);
    std::vector<AGD_TYPE *> data_ptr(n_patches, nullptr);
    struct BlendConfig blend_config = {};
    std::vector<struct BlendParams> blend_params(n_patches);
    AGD_TYPE ***blended;

    shottimes[0].push_back(-5);
    shottimes[0].push_back(-3);
    channels[0].resize(n_traces_per_patch[0]);
    std::fill(channels[0].begin(), channels[0].end(), 123);
    trace_types[0].push_back(AGDLive);
    trace_types[0].push_back(AGDLive);
    data[0].resize(n_traces_per_patch[0] * n_times_per_patch[0]);
    std::iota(data[0].begin(), data[0].end(), (AGD_TYPE)1.0);

    std::vector<AGD_TYPE> expected_chan0_count{
        weight[0],
        weight[1] + weight[0],
        weight[2] + weight[0] + weight[1],
        1 + weight[1] + weight[2],
        1 + weight[2] + 1,
        3,
        3,
        weight[2] + 2,
        weight[1] + 2,
        weight[0] + 2 * weight[2],
        2 * weight[1],
        2 * weight[0]};
    std::vector<AGD_TYPE> expected_chan0{
        1 * weight[0],
        2 * weight[1] + 61 * weight[0],
        3 * weight[2] + 11 * weight[0] + 62 * weight[1],
        4 + 12 * weight[1] + 63 * weight[2],
        5 + 13 * weight[2] + 64,
        6 + 14 + 65,
        7 + 15 + 66,
        8 * weight[2] + 16 + 67,
        9 * weight[1] + 17 + 68,
        10 * weight[0] + (18 + 69) * weight[2],
        (19 + 70) * weight[1],
        (20 + 71) * weight[0]};

    for (int i = 0; i < n_patches; i++) {
      shottimes_ptr[i] = &shottimes[i][0];
      channels_ptr[i] = &channels[i][0];
      trace_types_ptr[i] = &trace_types[i][0];
      data_ptr[i] = &data[i][0];
    }

    /* blend_config */
    EXPECT_EQ(set_blend_config(n_patches, n_traces_per_patch, n_times_per_patch,
                               shottimes_ptr.data(), channels_ptr.data(),
                               trace_types_ptr.data(), newcomm, &blend_config),
              0);
    EXPECT_EQ(blend_config.n_channels, 1);
    EXPECT_EQ(blend_config.n_mute_coords, 0);
    EXPECT_EQ(blend_config.unique_channel_values[0], 123);
    EXPECT_EQ(blend_config.channels_intervals[0].n_intervals, 1);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].start, -5);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].stop, 7);

    /* blend_params */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      EXPECT_EQ(set_blend_params(
                    n_traces_per_patch[patch_idx], n_times_per_patch[patch_idx],
                    shottimes[patch_idx].data(), channels[patch_idx].data(),
                    trace_types[patch_idx].data(), &blend_config,
                    blend_params.data() + patch_idx),
                0);
    }

    /* allocate blended */
    EXPECT_EQ(allocate_blended(&blend_config, &blended), 0);
    zero_blended(&blend_config, (AGD_TYPE *const *const *)blended);

    /* blend mean */
    blend_mean_forward(n_patches, data_ptr.data(), &blend_config,
                       blend_params.data(), taper_length,
                       (AGD_TYPE *const *const *)blended);

    for (int time_idx = 0; time_idx < 12; time_idx++) {
      EXPECT_FLOAT_EQ(
          blended[0][0][time_idx],
          expected_chan0[time_idx] / expected_chan0_count[time_idx]);
    }

    free_blended(&blend_config, &blended);
    for (auto &blend_param : blend_params) {
      free_blend_params(&blend_param);
    }
    free_blend_config(&blend_config);
  } else if (comm_rank == 1) {
    int n_patches = 1;
    int n_traces_per_patch[1] = {2};
    int n_times_per_patch[1] = {11};
    int taper_length = 3;
    AGD_TYPE weight[3] = {AGD_HALF * AGD_SQRT(AGD_TWO - AGD_SQRT(AGD_TWO)),
                          AGD_HALF * AGD_SQRT(AGD_TWO),
                          AGD_HALF * AGD_SQRT(AGD_TWO + AGD_SQRT(AGD_TWO))};
    std::vector<std::vector<long int>> shottimes(1);
    std::vector<long int *> shottimes_ptr(n_patches, nullptr);
    std::vector<std::vector<int>> channels(1);
    std::vector<int *> channels_ptr(n_patches, nullptr);
    std::vector<std::vector<enum AGDTraceType>> trace_types(1);
    std::vector<enum AGDTraceType *> trace_types_ptr(n_patches, nullptr);
    std::vector<std::vector<AGD_TYPE>> data(1);
    std::vector<AGD_TYPE *> data_ptr(n_patches, nullptr);
    struct BlendConfig blend_config = {};
    std::vector<struct BlendParams> blend_params(n_patches);
    AGD_TYPE ***blended;

    shottimes[0].push_back(-10);
    shottimes[0].push_back(-4);
    channels[0].resize(n_traces_per_patch[0]);
    std::fill(channels[0].begin(), channels[0].end(), 123);
    channels[0][0] = 456;
    trace_types[0].push_back(AGDLive);
    trace_types[0].push_back(AGDLive);
    data[0].resize(n_traces_per_patch[0] * n_times_per_patch[0]);
    std::iota(data[0].begin(), data[0].end(), (AGD_TYPE)50.0);

    std::vector<AGD_TYPE> expected_chan0_count{
        weight[1] + weight[0],
        weight[2] + weight[0] + weight[1],
        1 + weight[1] + weight[2],
        1 + weight[2] + 1,
        3,
        3,
        weight[2] + 2,
        weight[1] + 2,
        weight[0] + 2 * weight[2],
        2 * weight[1],
        2 * weight[0]};
    std::vector<AGD_TYPE> expected_chan0{
        2 * weight[1] + 61 * weight[0],
        3 * weight[2] + 11 * weight[0] + 62 * weight[1],
        4 + 12 * weight[1] + 63 * weight[2],
        5 + 13 * weight[2] + 64,
        6 + 14 + 65,
        7 + 15 + 66,
        8 * weight[2] + 16 + 67,
        9 * weight[1] + 17 + 68,
        10 * weight[0] + (18 + 69) * weight[2],
        (19 + 70) * weight[1],
        (20 + 71) * weight[0]};

    for (int i = 0; i < n_patches; i++) {
      shottimes_ptr[i] = &shottimes[i][0];
      channels_ptr[i] = &channels[i][0];
      trace_types_ptr[i] = &trace_types[i][0];
      data_ptr[i] = &data[i][0];
    }

    /* blend_config */
    EXPECT_EQ(set_blend_config(n_patches, n_traces_per_patch, n_times_per_patch,
                               shottimes_ptr.data(), channels_ptr.data(),
                               trace_types_ptr.data(), newcomm, &blend_config),
              0);
    EXPECT_EQ(blend_config.n_channels, 2);
    EXPECT_EQ(blend_config.n_mute_coords, 0);
    EXPECT_EQ(blend_config.unique_channel_values[0], 456);
    EXPECT_EQ(blend_config.unique_channel_values[1], 123);
    EXPECT_EQ(blend_config.channels_intervals[0].n_intervals, 1);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].start, -10);
    EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].stop, 1);
    EXPECT_EQ(blend_config.channels_intervals[1].n_intervals, 1);
    EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].start, -4);
    EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].stop, 7);

    /* blend_params */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      EXPECT_EQ(set_blend_params(
                    n_traces_per_patch[patch_idx], n_times_per_patch[patch_idx],
                    shottimes[patch_idx].data(), channels[patch_idx].data(),
                    trace_types[patch_idx].data(), &blend_config,
                    blend_params.data() + patch_idx),
                0);
    }

    /* allocate blended */
    EXPECT_EQ(allocate_blended(&blend_config, &blended), 0);
    zero_blended(&blend_config, (AGD_TYPE *const *const *)blended);

    /* blend mean */
    blend_mean_forward(n_patches, data_ptr.data(), &blend_config,
                       blend_params.data(), taper_length,
                       (AGD_TYPE *const *const *)blended);

    for (int time_idx = 0; time_idx < 11; time_idx++) {
      EXPECT_FLOAT_EQ(
          blended[1][0][time_idx],
          expected_chan0[time_idx] / expected_chan0_count[time_idx]);
    }
    for (int time_idx = 0; time_idx < 11; time_idx++) {
      EXPECT_FLOAT_EQ(blended[0][0][time_idx], data[0][time_idx]);
    }

    free_blended(&blend_config, &blended);
    for (auto &blend_param : blend_params) {
      free_blend_params(&blend_param);
    }
    free_blend_config(&blend_config);
  }

  MPI_Comm_free(&newcomm);
  // MPI_Finalize();
}

TEST(AGDBlend, TwoVolumesOneChannelSameOutMPI) {
  MPI_Comm newcomm;
  int comm_rank;
  int comm_size;

  // MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if (comm_size < 2) {
    fprintf(
        stderr,
        "WARNING: More than 2 ranks are required for this test. SKIPPING!\n");
    return;
  }
  MPI_Comm_split(MPI_COMM_WORLD, comm_rank < 2 ? 0 : MPI_UNDEFINED, 0,
                 &newcomm);
  MPI_Comm_rank(newcomm, &comm_rank);
  MPI_Comm_size(newcomm, &comm_size);

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
  int n_times = 3;

  shottimes = (long int **)malloc(1 * sizeof(long int *));
  channels = (int **)malloc(1 * sizeof(int *));
  trace_types = (AGDTraceType **)malloc(1 * sizeof(AGDTraceType *));
  data = (AGD_TYPE **)malloc(1 * sizeof(AGD_TYPE *));
  data_out = (AGD_TYPE **)malloc(1 * sizeof(AGD_TYPE *));
  shottimes[0] = (long int *)malloc(sizeof(long int));
  channels[0] = (int *)malloc(sizeof(int));
  trace_types[0] = (AGDTraceType *)malloc(sizeof(AGDTraceType));
  data[0] = (AGD_TYPE *)malloc(n_times * sizeof(AGD_TYPE));
  data_out[0] = (AGD_TYPE *)malloc(n_times * sizeof(AGD_TYPE));

  if (comm_rank == 0) {
    shottimes[0][0] = 25;
    channels[0][0] = 123;
    trace_types[0][0] = AGDLive;
    data[0][0] = 1;
    data[0][1] = 2;
    data[0][2] = 3;
  } else if (comm_rank == 1) {
    shottimes[0][0] = 27;
    channels[0][0] = 123;
    trace_types[0][0] = AGDLive;
    data[0][0] = 4;
    data[0][1] = 5;
    data[0][2] = 6;
  }

  shottimes_out = shottimes;
  channels_out = channels;
  trace_types_out = trace_types;

  if (comm_rank < 2) {
    EXPECT_EQ(agd_blend(n_patches, n_traces_per_patch, n_times_per_patch,
                        shottimes, channels, trace_types, data, blend_mode,
                        taper_length, n_patches_out, n_traces_per_patch_out,
                        n_times_per_patch_out, shottimes_out, channels_out,
                        trace_types_out, newcomm, data_out),
              0);
  }

  if (comm_rank == 0) {
    EXPECT_EQ(data_out[0][0], 1);
    EXPECT_EQ(data_out[0][1], 2);
    EXPECT_EQ(data_out[0][2], 7);
  } else if (comm_rank == 1) {
    EXPECT_EQ(data_out[0][0], 7);
    EXPECT_EQ(data_out[0][1], 5);
    EXPECT_EQ(data_out[0][2], 6);
  }

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

  MPI_Comm_free(&newcomm);
  // MPI_Finalize();
}

TEST(AGDeblend, Deblend) {
  MPI_Comm newcomm;
  int comm_rank;
  int comm_size;

  // MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if (comm_size < 2) {
    fprintf(
        stderr,
        "WARNING: More than 2 ranks are required for this test. SKIPPING!\n");
    return;
  }
  MPI_Comm_split(MPI_COMM_WORLD, comm_rank < 2 ? 0 : MPI_UNDEFINED, 0,
                 &newcomm);
  MPI_Comm_rank(newcomm, &comm_rank);
  MPI_Comm_size(newcomm, &comm_size);

  if (comm_rank < 2) {
    int n_patches = 1;
    int n_traces = 8;
    int n_times = 64;
    std::vector<int> n_dims = {2};
    std::vector<int> volumes = {0};
    std::vector<std::vector<int>> data_shapes = {{n_traces, n_times}};
    std::vector<int *> data_shapes_ptr = {&data_shapes[0][0]};
    std::vector<std::vector<int>> window_shapes = {{8, 32}};
    std::vector<int *> window_shapes_ptr = {&window_shapes[0][0]};
    std::vector<std::vector<int>> coords = {std::vector<int>(1)};
    std::vector<int *> coords_ptr = {&coords[0][0]};
    std::vector<std::vector<long int>> shottimes = {
        std::vector<long int>(n_traces)};
    std::vector<long int *> shottimes_ptr = {&shottimes[0][0]};
    std::vector<std::vector<int>> channels = {std::vector<int>(n_traces, 456)};
    std::vector<int *> channels_ptr = {&channels[0][0]};
    std::vector<std::vector<enum AGDTraceType>> trace_types_blend = {
        std::vector<enum AGDTraceType>(n_traces, AGDLive)};
    std::vector<enum AGDTraceType *> trace_types_blend_ptr = {
        &trace_types_blend[0][0]};
    std::vector<std::vector<enum AGDTraceType>> trace_types = {
        std::vector<enum AGDTraceType>(n_traces, AGDLive)};
    std::vector<enum AGDTraceType *> trace_types_ptr = {&trace_types[0][0]};
    std::vector<std::vector<AGD_TYPE>> true_data = {
        std::vector<AGD_TYPE>(12 * n_times)};
    std::vector<AGD_TYPE *> true_data_ptr = {&true_data[0][0]};
    std::vector<std::vector<AGD_TYPE>> joined_data = {
        std::vector<AGD_TYPE>(12 * n_times)};
    std::vector<AGD_TYPE *> joined_data_ptr = {&joined_data[0][0]};
    std::vector<std::vector<AGD_TYPE>> data = {
        std::vector<AGD_TYPE>(n_traces * n_times)};
    std::vector<AGD_TYPE *> data_ptr = {&data[0][0]};
    int *wavelet_lengths = nullptr;
    int **wavelet_idxs = nullptr;
    AGD_TYPE **wavelets = nullptr;

    if (comm_rank == 0) {
      coords[0][0] = 0;
      std::vector<long int> shottimes0 = {-1, 31, 57, 98, 118, 152, 179, 228};
      for (int trace_idx = 0; trace_idx < n_traces; trace_idx++) {
        shottimes[0][trace_idx] = shottimes0[trace_idx];
      }
      trace_types_blend[0][6] = AGDMissing;
      trace_types_blend[0][7] = AGDMissing;
    } else if (comm_rank == 1) {
      coords[0][0] = 1;
      trace_types_blend[0][0] = AGDMissing;
      trace_types_blend[0][1] = AGDMissing;
      std::vector<long int> shottimes1 = {118, 152, 179, 228,
                                          269, 292, 311, 365};
      for (int trace_idx = 0; trace_idx < n_traces; trace_idx++) {
        shottimes[0][trace_idx] = shottimes1[trace_idx];
      }
    }

    /* Make true data */
    std::vector<AGD_TYPE> time_values = {
        0.00239766,  -0.32091303, -0.02946537, -0.35292345, -0.04336042,
        -0.00480379, -0.19307142, -0.27653807, -0.05905185, 0.04624617,
        0.2010874,   0.41034904,  -0.2747792,  -0.09522921, 0.29506742,
        0.33675247,  0.03078658,  0.14770945,  0.21812879,  -0.37182982,
        -0.43089135, -0.41063197, 0.0058671,   -0.05369114, 0.38726809,
        0.2263909,   -0.11223266, 0.09796578,  -0.12503645, 0.21239305,
        0.0704477,   0.44837965,  -0.38498882, -0.42468094, 0.04627415,
        -0.23315139, 0.03198743,  -0.32517176, 0.28850815,  -0.18743274,
        0.31115991,  0.38803847,  0.07861607,  0.33403261,  -0.24939925,
        0.3799239,   -0.45458053, -0.26089656, -0.24444692, -0.21011644,
        0.32014846,  -0.06315092, 0.32757867,  -0.26892939, 0.17115324,
        -0.26522846, -0.60002929, 0.44351847,  0.37738665,  0.02023595,
        -0.23541187, 0.48802938,  0.40298129,  0.03456871};
    for (size_t trace_idx = 0; trace_idx < 12; trace_idx++) {
      for (size_t time_idx = 0; time_idx < time_values.size(); time_idx++) {
        true_data[0][trace_idx * n_times + time_idx] = time_values[time_idx];
      }
    }

    /* Extract this rank's section of true data */
    memcpy(&data[0][0], &true_data[0][comm_rank * 4 * n_times],
           (size_t)(n_traces * n_times) * sizeof(AGD_TYPE));

    /* Blend */
    EXPECT_EQ(agd_blend(n_patches, &n_traces, &n_times, shottimes_ptr.data(),
                        channels_ptr.data(), trace_types_blend_ptr.data(),
                        data_ptr.data(), AGDBlendSum, 0, n_patches, &n_traces,
                        &n_times, shottimes_ptr.data(), channels_ptr.data(),
                        trace_types_ptr.data(), newcomm, data_ptr.data()),
              0);

    /* Deblend */
    AGD_TYPE initial_factor = AGD_ONE;
    int n_its = 2000;
    EXPECT_EQ(agd_deblend(n_patches, volumes.data(), n_dims.data(),
                          window_shapes_ptr.data(), coords_ptr.data(),
                          data_shapes_ptr.data(), shottimes_ptr.data(),
                          channels_ptr.data(), trace_types_ptr.data(),
                          wavelet_lengths, wavelet_idxs, wavelets,
                          initial_factor, n_its, -1, newcomm, data_ptr.data()),
              0);

    /* Sum the different ranks' patches */
    memcpy(&joined_data[0][comm_rank * 4 * n_times], &data[0][0],
           (size_t)(n_traces * n_times) * sizeof(AGD_TYPE));
    MPI_Allreduce(MPI_IN_PLACE, &joined_data[0][0], 12 * n_times, AGD_MPI_TYPE,
                  MPI_SUM, newcomm);

    AGD_TYPE error_sq = AGD_ZERO;
    AGD_TYPE true_sq = AGD_ZERO;
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      for (int idx = 0; idx < 12 * n_times; idx++) {
        AGD_TYPE error =
            joined_data[patch_idx][idx] - true_data[patch_idx][idx];
        error_sq += error * error;
        true_sq += true_data[patch_idx][idx] * true_data[patch_idx][idx];
      }
    }
    EXPECT_EQ(error_sq < true_sq / 50, 1);
  }

  MPI_Comm_free(&newcomm);
  MPI_Finalize();
}
