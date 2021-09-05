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

TEST(Model, SetSpaceWindows) {
  int n_patches = 5;
  int volumes[5] = {0, 1, 2, 3, 3};
  int n_dims[4] = {2, 3, 3, 3};
  std::vector<std::vector<int>> coords = {{0}, {0, 0}, {0, 0}, {0, 0}, {1, 1}};
  std::vector<int *> coords_ptr = {&coords[0][0], &coords[1][0], &coords[2][0],
                                   &coords[3][0], &coords[4][0]};
  std::vector<std::vector<int>> data_shapes = {
      {4, 8}, {3, 5, 7}, {3, 7, 7}, {4, 2, 4}, {6, 2, 4}};
  std::vector<int *> data_shapes_ptr = {&data_shapes[0][0], &data_shapes[1][0],
                                        &data_shapes[2][0], &data_shapes[3][0],
                                        &data_shapes[4][0]};
  std::vector<std::vector<int>> window_shapes = {
      {2, 8}, {3, 2, 6}, {2, 4, 6}, {4, 2, 4}};
  std::vector<int *> window_shapes_ptr = {
      &window_shapes[0][0], &window_shapes[1][0], &window_shapes[2][0],
      &window_shapes[3][0]};
  int n_windows;
  struct ModelWindow *windows;

  EXPECT_EQ(set_space_windows(n_patches, volumes, n_dims,
                              window_shapes_ptr.data(), coords_ptr.data(),
                              data_shapes_ptr.data(), &n_windows, &windows),
            0);

  /* Volume 0:
   * 0 01 12 2
   *
   * Volume 1:
   * 3  3  3
   * 34 34 34
   * 45 45 45
   * 56 56 56
   * 6  6  6
   *
   * Volume 2:
   * 7  79   9
   * 7  79   9
   * 7  79   9
   * 78 789X 9X
   * 78 789X 9X
   * 8  8X   X
   * 8  8X   X
   *
   * Volume 3:
   * 11 11 11 11
   * 11 11 11 11
   *             12 12 12/13 12/13 13 13
   *             12 12 12/13 12/13 13 13
   */

  EXPECT_EQ(n_windows, 3 + 4 + 4 + 3);

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
  EXPECT_EQ(window.use_taper[1], 0);

  window = windows[3].space_window;
  EXPECT_EQ(window.patch_idx, 1);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 6);
  EXPECT_EQ(window.coords[0], 0);
  EXPECT_EQ(window.coords[1], 0);
  EXPECT_EQ(window.space_shape[0], 3);
  EXPECT_EQ(window.space_shape[1], 2);
  EXPECT_EQ(window.taper_length[1], 1);
  EXPECT_EQ(window.use_taper[0], 0);
  EXPECT_EQ(window.use_taper[1], 0);
  EXPECT_EQ(window.use_taper[2], 0);
  EXPECT_EQ(window.use_taper[3], 1);

  window = windows[4].space_window;
  EXPECT_EQ(window.patch_idx, 1);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 6);
  EXPECT_EQ(window.coords[0], 0);
  EXPECT_EQ(window.coords[1], 1);
  EXPECT_EQ(window.space_shape[0], 3);
  EXPECT_EQ(window.space_shape[1], 2);
  EXPECT_EQ(window.taper_length[1], 1);
  EXPECT_EQ(window.use_taper[0], 0);
  EXPECT_EQ(window.use_taper[1], 0);
  EXPECT_EQ(window.use_taper[2], 1);
  EXPECT_EQ(window.use_taper[3], 1);

  window = windows[5].space_window;
  EXPECT_EQ(window.patch_idx, 1);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 6);
  EXPECT_EQ(window.coords[0], 0);
  EXPECT_EQ(window.coords[1], 2);
  EXPECT_EQ(window.space_shape[0], 3);
  EXPECT_EQ(window.space_shape[1], 2);
  EXPECT_EQ(window.taper_length[1], 1);
  EXPECT_EQ(window.use_taper[0], 0);
  EXPECT_EQ(window.use_taper[1], 0);
  EXPECT_EQ(window.use_taper[2], 1);
  EXPECT_EQ(window.use_taper[3], 1);

  window = windows[6].space_window;
  EXPECT_EQ(window.patch_idx, 1);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 6);
  EXPECT_EQ(window.coords[0], 0);
  EXPECT_EQ(window.coords[1], 3);
  EXPECT_EQ(window.space_shape[0], 3);
  EXPECT_EQ(window.space_shape[1], 2);
  EXPECT_EQ(window.taper_length[1], 1);
  EXPECT_EQ(window.use_taper[0], 0);
  EXPECT_EQ(window.use_taper[1], 0);
  EXPECT_EQ(window.use_taper[2], 1);
  EXPECT_EQ(window.use_taper[3], 0);

  window = windows[7].space_window;
  EXPECT_EQ(window.patch_idx, 2);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 10);
  EXPECT_EQ(window.coords[0], 0);
  EXPECT_EQ(window.coords[1], 0);
  EXPECT_EQ(window.space_shape[0], 2);
  EXPECT_EQ(window.space_shape[1], 5);
  EXPECT_EQ(window.taper_length[0], 1);
  EXPECT_EQ(window.taper_length[1], 2);
  EXPECT_EQ(window.use_taper[0], 0);
  EXPECT_EQ(window.use_taper[1], 1);
  EXPECT_EQ(window.use_taper[2], 0);
  EXPECT_EQ(window.use_taper[3], 1);

  window = windows[8].space_window;
  EXPECT_EQ(window.patch_idx, 2);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 8);
  EXPECT_EQ(window.coords[0], 0);
  EXPECT_EQ(window.coords[1], 3);
  EXPECT_EQ(window.space_shape[0], 2);
  EXPECT_EQ(window.space_shape[1], 4);
  EXPECT_EQ(window.taper_length[0], 1);
  EXPECT_EQ(window.taper_length[1], 2);
  EXPECT_EQ(window.use_taper[0], 0);
  EXPECT_EQ(window.use_taper[1], 1);
  EXPECT_EQ(window.use_taper[2], 1);
  EXPECT_EQ(window.use_taper[3], 0);

  window = windows[9].space_window;
  EXPECT_EQ(window.patch_idx, 2);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 10);
  EXPECT_EQ(window.coords[0], 1);
  EXPECT_EQ(window.coords[1], 0);
  EXPECT_EQ(window.space_shape[0], 2);
  EXPECT_EQ(window.space_shape[1], 5);
  EXPECT_EQ(window.taper_length[0], 1);
  EXPECT_EQ(window.taper_length[1], 2);
  EXPECT_EQ(window.use_taper[0], 1);
  EXPECT_EQ(window.use_taper[1], 0);
  EXPECT_EQ(window.use_taper[2], 0);
  EXPECT_EQ(window.use_taper[3], 1);

  window = windows[10].space_window;
  EXPECT_EQ(window.patch_idx, 2);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 8);
  EXPECT_EQ(window.coords[0], 1);
  EXPECT_EQ(window.coords[1], 3);
  EXPECT_EQ(window.space_shape[0], 2);
  EXPECT_EQ(window.space_shape[1], 4);
  EXPECT_EQ(window.taper_length[0], 1);
  EXPECT_EQ(window.taper_length[1], 2);
  EXPECT_EQ(window.use_taper[0], 1);
  EXPECT_EQ(window.use_taper[1], 0);
  EXPECT_EQ(window.use_taper[2], 1);
  EXPECT_EQ(window.use_taper[3], 0);

  window = windows[11].space_window;
  EXPECT_EQ(window.patch_idx, 3);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 8);
  EXPECT_EQ(window.coords[0], 0);
  EXPECT_EQ(window.coords[1], 0);
  EXPECT_EQ(window.space_shape[0], 4);
  EXPECT_EQ(window.space_shape[1], 2);
  EXPECT_EQ(window.taper_length[0], 2);
  EXPECT_EQ(window.taper_length[1], 1);
  EXPECT_EQ(window.use_taper[0], 0);
  EXPECT_EQ(window.use_taper[1], 1);
  EXPECT_EQ(window.use_taper[2], 0);
  EXPECT_EQ(window.use_taper[3], 1);

  window = windows[12].space_window;
  EXPECT_EQ(window.patch_idx, 4);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 8);
  EXPECT_EQ(window.coords[0], 0);
  EXPECT_EQ(window.coords[1], 0);
  EXPECT_EQ(window.space_shape[0], 4);
  EXPECT_EQ(window.space_shape[1], 2);
  EXPECT_EQ(window.taper_length[0], 2);
  EXPECT_EQ(window.taper_length[1], 1);
  EXPECT_EQ(window.use_taper[0], 1);
  EXPECT_EQ(window.use_taper[1], 1);
  EXPECT_EQ(window.use_taper[2], 1);
  EXPECT_EQ(window.use_taper[3], 0);

  window = windows[13].space_window;
  EXPECT_EQ(window.patch_idx, 4);
  EXPECT_EQ(window.n_dims, 3);
  EXPECT_EQ(window.n_traces, 8);
  EXPECT_EQ(window.coords[0], 2);
  EXPECT_EQ(window.coords[1], 0);
  EXPECT_EQ(window.space_shape[0], 4);
  EXPECT_EQ(window.space_shape[1], 2);
  EXPECT_EQ(window.taper_length[0], 2);
  EXPECT_EQ(window.taper_length[1], 1);
  EXPECT_EQ(window.use_taper[0], 1);
  EXPECT_EQ(window.use_taper[1], 0);
  EXPECT_EQ(window.use_taper[2], 1);
  EXPECT_EQ(window.use_taper[3], 0);

  for (int window_idx = 0; window_idx < n_windows; window_idx++) {
    struct ModelWindow *const window = windows + window_idx;
    free_space_window(&(window->space_window));
  }
  free(windows);
}

TEST(Model, SetSpaceWindowsLarger) {
  /* Check that extra traces are correctly added to windows when the window
   * length does not divide evenly into the data length.
   * With a window shape of 32, the taper length is 16. If the data length
   * is 78, we can fit 3 windows, but will be 14 traces short of covering
   * the data length. We thus need a base window length of 36
   * (16 left taper, 4 middle, 16 right taper), and the first 2 windows
   * will need to be one larger (37). */
  int n_patches = 1;
  int volumes[1] = {0};
  int n_dims[1] = {2};
  std::vector<std::vector<int>> coords = {{0}};
  std::vector<int *> coords_ptr = {&coords[0][0]};
  std::vector<std::vector<int>> data_shapes = {{78, 8}};
  std::vector<int *> data_shapes_ptr = {&data_shapes[0][0]};
  std::vector<std::vector<int>> window_shapes = {{32, 8}};
  std::vector<int *> window_shapes_ptr = {&window_shapes[0][0]};
  int n_windows;
  struct ModelWindow *windows;

  EXPECT_EQ(set_space_windows(n_patches, volumes, n_dims,
                              window_shapes_ptr.data(), coords_ptr.data(),
                              data_shapes_ptr.data(), &n_windows, &windows),
            0);

  EXPECT_EQ(n_windows, 3);

  struct SpaceWindow window;

  window = windows[0].space_window;
  EXPECT_EQ(window.patch_idx, 0);
  EXPECT_EQ(window.n_dims, 2);
  EXPECT_EQ(window.n_traces, 37);
  EXPECT_EQ(window.coords[0], 0);
  EXPECT_EQ(window.space_shape[0], 37);
  EXPECT_EQ(window.taper_length[0], 16);
  EXPECT_EQ(window.use_taper[0], 0);
  EXPECT_EQ(window.use_taper[1], 1);

  window = windows[1].space_window;
  EXPECT_EQ(window.patch_idx, 0);
  EXPECT_EQ(window.n_dims, 2);
  EXPECT_EQ(window.n_traces, 37);
  EXPECT_EQ(window.coords[0], 21);
  EXPECT_EQ(window.space_shape[0], 37);
  EXPECT_EQ(window.taper_length[0], 16);
  EXPECT_EQ(window.use_taper[0], 1);
  EXPECT_EQ(window.use_taper[1], 1);

  window = windows[2].space_window;
  EXPECT_EQ(window.patch_idx, 0);
  EXPECT_EQ(window.n_dims, 2);
  EXPECT_EQ(window.n_traces, 36);
  EXPECT_EQ(window.coords[0], 42);
  EXPECT_EQ(window.space_shape[0], 36);
  EXPECT_EQ(window.taper_length[0], 16);
  EXPECT_EQ(window.use_taper[0], 1);
  EXPECT_EQ(window.use_taper[1], 0);

  for (int window_idx = 0; window_idx < n_windows; window_idx++) {
    struct ModelWindow *const window = windows + window_idx;
    free_space_window(&(window->space_window));
  }
  free(windows);
}

TEST(Model, GetWindowArrays) {
  int n_patches = 1;
  int volumes[1] = {0};
  int n_dims[1] = {3};
  std::vector<std::vector<int>> coords = {{0, 0}};
  std::vector<int *> coords_ptr = {&coords[0][0]};
  std::vector<std::vector<int>> data_shapes = {{3, 7, 7}};
  std::vector<int *> data_shapes_ptr = {&data_shapes[0][0]};
  std::vector<std::vector<int>> window_shapes = {{2, 4, 6}};
  std::vector<int *> window_shapes_ptr = {&window_shapes[0][0]};
  std::vector<std::vector<long int>> shottimes({std::vector<long int>(3 * 7)});
  std::vector<long int *> shottimes_ptr = {&shottimes[0][0]};
  std::vector<std::vector<int>> channels({std::vector<int>(3 * 7)});
  std::vector<int *> channels_ptr = {&channels[0][0]};
  std::vector<std::vector<enum AGDTraceType>> trace_types(
      {std::vector<enum AGDTraceType>(3 * 7)});
  std::vector<enum AGDTraceType *> trace_types_ptr = {&trace_types[0][0]};
  int n_windows;
  struct ModelWindow *windows;

  /*
   0  1  2  3  4  5  6
   7  8  9 10 11 12 13
  14 15 16 17 18 19 20
  */

  srand(1);
  std::iota(shottimes[0].begin(), shottimes[0].end(), 0);
  std::fill(channels[0].begin(), channels[0].end(), 0);
  std::fill(trace_types[0].begin(), trace_types[0].end(), AGDLive);

  EXPECT_EQ(set_space_windows(n_patches, volumes, n_dims,
                              window_shapes_ptr.data(), coords_ptr.data(),
                              data_shapes_ptr.data(), &n_windows, &windows),
            0);
  EXPECT_EQ(n_windows, 4);

  for (int window_idx = 0; window_idx < n_windows; window_idx++) {
    struct WindowArrays window_arrays = ZERO_INIT;
    struct ModelWindow *const window = windows + window_idx;
    EXPECT_EQ(
        set_window_arrays(&(window->space_window), data_shapes_ptr.data(),
                          shottimes_ptr.data(), channels_ptr.data(),
                          trace_types_ptr.data(), nullptr, &window_arrays),
        0);

    if (window_idx == 0) {
      EXPECT_EQ(window_arrays.shottimes[0], 0);
      EXPECT_EQ(window_arrays.shottimes[1], 1);
      EXPECT_EQ(window_arrays.shottimes[2], 2);
      EXPECT_EQ(window_arrays.shottimes[3], 3);
      EXPECT_EQ(window_arrays.shottimes[4], 4);
      EXPECT_EQ(window_arrays.shottimes[5], 7);
      EXPECT_EQ(window_arrays.shottimes[6], 8);
      EXPECT_EQ(window_arrays.shottimes[7], 9);
      EXPECT_EQ(window_arrays.shottimes[8], 10);
      EXPECT_EQ(window_arrays.shottimes[9], 11);
    } else if (window_idx == 1) {
      EXPECT_EQ(window_arrays.shottimes[0], 3);
      EXPECT_EQ(window_arrays.shottimes[1], 4);
      EXPECT_EQ(window_arrays.shottimes[2], 5);
      EXPECT_EQ(window_arrays.shottimes[3], 6);
      EXPECT_EQ(window_arrays.shottimes[4], 10);
      EXPECT_EQ(window_arrays.shottimes[5], 11);
      EXPECT_EQ(window_arrays.shottimes[6], 12);
      EXPECT_EQ(window_arrays.shottimes[7], 13);
    } else if (window_idx == 2) {
      EXPECT_EQ(window_arrays.shottimes[0], 7);
      EXPECT_EQ(window_arrays.shottimes[1], 8);
      EXPECT_EQ(window_arrays.shottimes[2], 9);
      EXPECT_EQ(window_arrays.shottimes[3], 10);
      EXPECT_EQ(window_arrays.shottimes[4], 11);
      EXPECT_EQ(window_arrays.shottimes[5], 14);
      EXPECT_EQ(window_arrays.shottimes[6], 15);
      EXPECT_EQ(window_arrays.shottimes[7], 16);
      EXPECT_EQ(window_arrays.shottimes[8], 17);
      EXPECT_EQ(window_arrays.shottimes[9], 18);
    } else {
      EXPECT_EQ(window_arrays.shottimes[0], 10);
      EXPECT_EQ(window_arrays.shottimes[1], 11);
      EXPECT_EQ(window_arrays.shottimes[2], 12);
      EXPECT_EQ(window_arrays.shottimes[3], 13);
      EXPECT_EQ(window_arrays.shottimes[4], 17);
      EXPECT_EQ(window_arrays.shottimes[5], 18);
      EXPECT_EQ(window_arrays.shottimes[6], 19);
      EXPECT_EQ(window_arrays.shottimes[7], 20);
    }

    free_window_arrays(&window_arrays);
  }

  for (int window_idx = 0; window_idx < n_windows; window_idx++) {
    struct ModelWindow *const window = windows + window_idx;
    free_space_window(&(window->space_window));
  }
  free(windows);
}

TEST(Model, GetStepsize) {
  int n_patches = 1;
  int volumes[1] = {0};
  std::vector<int> n_dims = {2};
  std::vector<std::vector<int>> data_shapes = {{16, 16}};
  std::vector<int *> data_shapes_ptr = {&data_shapes[0][0]};
  std::vector<std::vector<int>> window_shapes = {{4, 8}};
  std::vector<int *> window_shapes_ptr = {&window_shapes[0][0]};
  std::vector<std::vector<int>> coords = {{0}};
  std::vector<int *> coords_ptr = {&coords[0][0]};
  std::vector<std::vector<long int>> shottimes = {std::vector<long int>(16)};
  std::vector<long int *> shottimes_ptr = {&shottimes[0][0]};
  std::vector<std::vector<int>> channels = {std::vector<int>(16)};
  std::vector<int *> channels_ptr = {&channels[0][0]};
  std::vector<std::vector<enum AGDTraceType>> trace_types = {
      std::vector<enum AGDTraceType>(16, AGDLive)};
  std::vector<enum AGDTraceType *> trace_types_ptr = {&trace_types[0][0]};
  struct Model model = {};
  AGD_TYPE ***blended;

  for (int i = 0; i < 16; i++) {
    shottimes[0][i] = (16 + 8) * i;
  }

  EXPECT_EQ(
      set_model(n_patches, volumes, n_dims.data(), window_shapes_ptr.data(),
                coords_ptr.data(), data_shapes_ptr.data(), shottimes_ptr.data(),
                channels_ptr.data(), trace_types_ptr.data(), nullptr, nullptr,
                nullptr, &model),
      0);

  EXPECT_EQ(allocate_blended(&(model.blend_config), &blended), 0);
  EXPECT_EQ(get_stepsize(&model, blended), AGD_HALF);
  for (int time_idx = 0; time_idx < 16 * (16 + 8); time_idx++) {
    EXPECT_FLOAT_EQ(blended[0][0][time_idx], AGD_ONE);
  }

  free_blended(&(model.blend_config), &blended);
  free_model(&model);

  for (int i = 0; i < 16; i++) {
    shottimes[0][i] = (16 + 7) * i;
  }

  EXPECT_EQ(
      set_model(n_patches, volumes, n_dims.data(), window_shapes_ptr.data(),
                coords_ptr.data(), data_shapes_ptr.data(), shottimes_ptr.data(),
                channels_ptr.data(), trace_types_ptr.data(), nullptr, nullptr,
                nullptr, &model),
      0);

  EXPECT_EQ(allocate_blended(&(model.blend_config), &blended), 0);
  EXPECT_EQ(get_stepsize(&model, blended), AGD_HALF / AGD_TWO);

  free_blended(&(model.blend_config), &blended);
  free_model(&model);
}

TEST(Model, Forward) {
  int n_patches = 4;
  int volumes[4] = {0, 1, 2, 0};
  std::vector<int> n_dims = {2, 3, 2, 2};
  std::vector<std::vector<int>> data_shapes = {
      {4, 4}, {2, 3, 5}, {4, 4}, {4, 4}};
  std::vector<int *> data_shapes_ptr = {&data_shapes[0][0], &data_shapes[1][0],
                                        &data_shapes[2][0], &data_shapes[3][0]};
  /* patch 0+3 will have three windows, with 8 time samples each
   * patch 1 will have two windows, with 10 time samples each */
  std::vector<std::vector<int>> window_shapes = {{2, 4}, {2, 2, 4}, {4, 4}};
  std::vector<int *> window_shapes_ptr = {
      &window_shapes[0][0], &window_shapes[1][0], &window_shapes[2][0]};
  std::vector<std::vector<int>> coords = {{0}, {0, 0}, {0}, {1}};
  std::vector<int *> coords_ptr = {&coords[0][0], &coords[1][0], &coords[2][0],
                                   &coords[3][0]};
  std::vector<std::vector<long int>> shottimes = {{1, 15, 51, 21},
                                                  {7, 27, 35, 43, 100, 57},
                                                  {1, 2, 3, 4},
                                                  {21, 65, 71, 77}};
  std::vector<long int *> shottimes_ptr = {&shottimes[0][0], &shottimes[1][0],
                                           &shottimes[2][0], &shottimes[3][0]};
  std::vector<std::vector<int>> channels = {
      {1, 1, 1, 1}, {1, 1, 1, 1, 2, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
  std::vector<int *> channels_ptr = {&channels[0][0], &channels[1][0],
                                     &channels[2][0], &channels[3][0]};
  std::vector<std::vector<enum AGDTraceType>> trace_types = {
      std::vector<enum AGDTraceType>(4, AGDLive),
      std::vector<enum AGDTraceType>(6, AGDLive),
      std::vector<enum AGDTraceType>(4, AGDMissing),
      std::vector<enum AGDTraceType>(4, AGDLive)};
  std::vector<enum AGDTraceType *> trace_types_ptr = {
      &trace_types[0][0], &trace_types[1][0], &trace_types[2][0],
      &trace_types[3][0]};
  int wavelet_lengths[1] = {1};
  std::vector<std::vector<int>> wavelet_idxs = {
      {0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  std::vector<int *> wavelet_idxs_ptr = {
      &wavelet_idxs[0][0], &wavelet_idxs[1][0], &wavelet_idxs[2][0],
      &wavelet_idxs[3][0]};
  std::vector<std::vector<AGD_TYPE>> wavelets = {{AGD_ONE}};
  std::vector<AGD_TYPE *> wavelets_ptr = {&wavelets[0][0]};
  struct Model model = {};
  AGD_FFTW_COMPLEX **model_fk;
  AGD_TYPE ***blended;

  EXPECT_EQ(
      set_model(n_patches, volumes, n_dims.data(), window_shapes_ptr.data(),
                coords_ptr.data(), data_shapes_ptr.data(), shottimes_ptr.data(),
                channels_ptr.data(), trace_types_ptr.data(), wavelet_lengths,
                wavelet_idxs_ptr.data(), wavelets_ptr.data(), &model),
      0);
  EXPECT_EQ(model.n_windows, 8);
  EXPECT_EQ(model.windows[0].space_window.patch_idx, 0);
  EXPECT_EQ(model.windows[1].space_window.patch_idx, 0);
  EXPECT_EQ(model.windows[2].space_window.patch_idx, 0);
  EXPECT_EQ(model.windows[3].space_window.patch_idx, 1);
  EXPECT_EQ(model.windows[4].space_window.patch_idx, 1);
  EXPECT_EQ(model.windows[5].space_window.patch_idx, 3);
  EXPECT_EQ(model.windows[6].space_window.patch_idx, 3);
  EXPECT_EQ(model.windows[7].space_window.patch_idx, 3);
  EXPECT_EQ(model.blend_config.n_channels, 2);
  EXPECT_EQ(model.blend_config.channels_intervals[0].n_intervals, 1);
  EXPECT_EQ(model.blend_config.channels_intervals[0].intervals[0].start, -1);
  EXPECT_EQ(model.blend_config.channels_intervals[0].intervals[0].stop, 83);
  EXPECT_EQ(model.blend_config.channels_intervals[1].n_intervals, 1);
  EXPECT_EQ(model.blend_config.channels_intervals[1].intervals[0].start, 98);
  EXPECT_EQ(model.blend_config.channels_intervals[1].intervals[0].stop, 108);

  EXPECT_EQ(allocate_model_fk(&model, &model_fk), 0);
  for (int window_idx = 0; window_idx < model.n_windows; window_idx++) {
    struct ModelWindow const *const window = &(model.windows[window_idx]);
    struct FKWindowConfig const *const fk_config =
        model.fk_configs.window_configs + window->fk_config_idx;
    for (int time_window_idx = 0; time_window_idx < fk_config->n_time_windows;
         time_window_idx++) {
      model_fk[window_idx]
              [time_window_idx * fk_config->n_elements_padded_window_fk][0] =
                  AGD_SQRT(fk_config->n_elements_padded_window);
    }
  }
  EXPECT_EQ(allocate_blended(&(model.blend_config), &blended), 0);
  zero_blended(&(model.blend_config), blended);
  model_forward(model_fk, &model, blended);

  EXPECT_NEAR(blended[0][0][0], AGD_ZERO, 1e-15);
  EXPECT_FLOAT_EQ(blended[0][0][1], AGD_HALF);
  for (int time_idx = 2; time_idx < 83; time_idx++) {
    EXPECT_FLOAT_EQ(blended[0][0][time_idx], AGD_ONE);
  }
  EXPECT_FLOAT_EQ(blended[0][0][83], AGD_HALF);

  EXPECT_NEAR(blended[1][0][0], AGD_ZERO, 1e-15);
  EXPECT_FLOAT_EQ(blended[1][0][1], AGD_HALF);
  for (int time_idx = 2; time_idx < 9; time_idx++) {
    EXPECT_FLOAT_EQ(blended[1][0][time_idx], AGD_ONE);
  }
  EXPECT_FLOAT_EQ(blended[1][0][9], AGD_HALF);

  free_blended(&(model.blend_config), &blended);
  free_model_fk(model.n_windows, &model_fk);
  free_model(&model);
}

TEST(Model, DotTest) {
  srand(1);
  for (int rand_repeat = 0; rand_repeat < 5; rand_repeat++) {
    int n_patches = 1 + rand() % 3;
    std::vector<int> volumes(n_patches);
    std::vector<int> n_dims(n_patches);
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
    std::vector<int> wavelet_lengths;
    std::vector<std::vector<int>> wavelet_idxs(n_patches);
    std::vector<int *> wavelet_idxs_ptr(n_patches);
    std::vector<std::vector<AGD_TYPE>> wavelets;
    std::vector<AGD_TYPE *> wavelets_ptr;
    struct Model model = ZERO_INIT;
    AGD_FFTW_COMPLEX **x;
    AGD_FFTW_COMPLEX **xp;
    AGD_TYPE ***y;
    AGD_TYPE ***yp;
    std::vector<int> n_traces(n_patches);
    int n_times = 0;

    /* Set data_shapes, window_shapes, n_traces, and n_times */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      data_shapes[patch_idx].resize(n_dims[patch_idx]);
      window_shapes[patch_idx].resize(n_dims[patch_idx]);
      coords[patch_idx].resize(n_dims[patch_idx] - 1);
      int n_times_preconv = 2 + rand() % 7;
      int n_times_conv = n_times_preconv + rand() % 4;
      for (int dim_idx = 0; dim_idx < n_dims[patch_idx]; dim_idx++) {
        if (dim_idx == n_dims[patch_idx] - 1) {
          data_shapes[patch_idx][dim_idx] = n_times_conv;
        } else {
          data_shapes[patch_idx][dim_idx] = 2 + rand() % 7;
        }
        window_shapes[patch_idx][dim_idx] =
            2 + rand() % (data_shapes[patch_idx][dim_idx] - 1);
        /* round down to even if not equal to data shape or time dim */
        if (window_shapes[patch_idx][dim_idx] !=
                data_shapes[patch_idx][dim_idx] ||
            dim_idx == n_dims[patch_idx] - 1) {
          window_shapes[patch_idx][dim_idx] =
              (window_shapes[patch_idx][dim_idx] / 2) * 2;
        }
      }
      n_traces[patch_idx] = 1;
      for (int dim_idx = 0; dim_idx < n_dims[patch_idx] - 1; dim_idx++) {
        n_traces[patch_idx] *= data_shapes[patch_idx][dim_idx];
      }
      n_times +=
          n_traces[patch_idx] * data_shapes[patch_idx][n_dims[patch_idx] - 1];

      int wavelet_length = n_times_conv - n_times_preconv + 1;
      /* Create some new wavelets of this length */
      int n_wavelets_to_create = 1 + rand() % (n_traces[patch_idx] - 1);
      for (int wavelet_idx = 0; wavelet_idx < n_wavelets_to_create;
           wavelet_idx++) {
        std::vector<AGD_TYPE> wavelet(wavelet_length);
        std::generate(wavelet.begin(), wavelet.end(), [] {
          return (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX - AGD_HALF;
        });
        wavelets.push_back(wavelet);
        wavelet_lengths.push_back(wavelet_length);
      }

      /* Make a list of idxs of wavelets of the right length */
      std::vector<int> possible_idxs;
      for (int wavelet_idx = 0;
           static_cast<size_t>(wavelet_idx) < wavelets.size(); wavelet_idx++) {
        if (wavelet_lengths[wavelet_idx] == wavelet_length) {
          possible_idxs.push_back(wavelet_idx);
        }
      }
      size_t n_possible = possible_idxs.size();

      /* Choose a wavelet from the possible idxs for each trace */
      for (int trace_idx = 0; trace_idx < n_traces[patch_idx]; trace_idx++) {
        wavelet_idxs[patch_idx].push_back(possible_idxs[rand() % n_possible]);
      }
    }

    /* Set shottimes, channels, trace_types */
    int shottime_max = n_times / 3;
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
      std::generate(trace_types[patch_idx].begin(),
                    trace_types[patch_idx].end(),
                    []() { return (enum AGDTraceType)(rand() % 3); });
    }

    /* Get pointers */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      data_shapes_ptr[patch_idx] = &data_shapes[patch_idx][0];
      window_shapes_ptr[patch_idx] = &window_shapes[patch_idx][0];
      coords_ptr[patch_idx] = &coords[patch_idx][0];
      shottimes_ptr[patch_idx] = &shottimes[patch_idx][0];
      channels_ptr[patch_idx] = &channels[patch_idx][0];
      trace_types_ptr[patch_idx] = &trace_types[patch_idx][0];
      wavelet_idxs_ptr[patch_idx] = &wavelet_idxs[patch_idx][0];
    }
    for (int wavelet_idx = 0;
         static_cast<size_t>(wavelet_idx) < wavelets.size(); wavelet_idx++) {
      wavelets_ptr.push_back(&wavelets[wavelet_idx][0]);
    }

    int use_wavelets = rand() > 0.5;
    int *wavelet_lengths_arg;
    int **wavelet_idxs_arg;
    AGD_TYPE **wavelets_arg;
    if (use_wavelets) {
      wavelet_lengths_arg = wavelet_lengths.data();
      wavelet_idxs_arg = wavelet_idxs_ptr.data();
      wavelets_arg = wavelets_ptr.data();
    } else {
      wavelet_lengths_arg = nullptr;
      wavelet_idxs_arg = nullptr;
      wavelets_arg = nullptr;
    }

    std::iota(volumes.begin(), volumes.end(), 0);

    /* Get model */
    EXPECT_EQ(
        set_model(n_patches, volumes.data(), n_dims.data(),
                  window_shapes_ptr.data(), coords_ptr.data(),
                  data_shapes_ptr.data(), shottimes_ptr.data(),
                  channels_ptr.data(), trace_types_ptr.data(),
                  wavelet_lengths_arg, wavelet_idxs_arg, wavelets_arg, &model),
        0);
    EXPECT_EQ(allocate_model_fk(&model, &x), 0);
    EXPECT_EQ(allocate_model_fk(&model, &xp), 0);
    EXPECT_EQ(allocate_blended(&(model.blend_config), &y), 0);
    EXPECT_EQ(allocate_blended(&(model.blend_config), &yp), 0);

    for (int window_idx = 0; window_idx < model.n_windows; window_idx++) {
      struct ModelWindow const *const window = model.windows + window_idx;
      struct FKWindowConfig const *const fk_config =
          model.fk_configs.window_configs + window->fk_config_idx;

      /* Make x by transforming random real numbers */
      for (int idx = 0; idx < fk_config->n_time_windows *
                                  fk_config->n_elements_padded_window;
           idx++) {
        static_cast<AGD_TYPE *>(model.temp_2[0])[idx] =
            (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX - AGD_HALF;
      }
      AGD_FFTW_EXECUTE(fk_config->adj_plan);
      memcpy(x[window_idx], model.temp_1[0],
             fk_config->n_time_windows *
                 fk_config->n_elements_padded_window_fk *
                 sizeof(AGD_FFTW_COMPLEX));

      memset(xp[window_idx], 0,
             fk_config->n_time_windows *
                 fk_config->n_elements_padded_window_fk *
                 sizeof(AGD_FFTW_COMPLEX));
    }

    for (int channel_idx = 0; channel_idx < model.blend_config.n_channels;
         channel_idx++) {
      struct ChannelIntervals const *const channel_intervals =
          model.blend_config.channels_intervals + channel_idx;
      for (int interval_idx = 0; interval_idx < channel_intervals->n_intervals;
           interval_idx++) {
        struct Interval const *const interval =
            channel_intervals->intervals + interval_idx;
        int n_times = interval->stop - interval->start;
        for (int time_idx = 0; time_idx < n_times; time_idx++) {
          yp[channel_idx][interval_idx][time_idx] =
              (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX;
        }
      }
    }

    apply_mute(&(model.blend_config), yp);
    zero_blended(&(model.blend_config), y);

    model_forward(x, &model, y);
    model_adjoint_update(yp, &model, -AGD_ONE, xp);

    AGD_TYPE sum_x = AGD_ZERO;
    for (int window_idx = 0; window_idx < model.n_windows; window_idx++) {
      struct ModelWindow const *const window = model.windows + window_idx;
      struct FKWindowConfig const *const fk_config =
          model.fk_configs.window_configs + window->fk_config_idx;
      for (int idx = 0; idx < fk_config->n_time_windows *
                                  fk_config->n_elements_padded_window_fk;
           idx++) {
        sum_x += x[window_idx][idx][0] * xp[window_idx][idx][0] +
                 x[window_idx][idx][1] * xp[window_idx][idx][1];
      }
    }

    AGD_TYPE sum_y = AGD_ZERO;
    for (int channel_idx = 0; channel_idx < model.blend_config.n_channels;
         channel_idx++) {
      struct ChannelIntervals const *const channel_intervals =
          model.blend_config.channels_intervals + channel_idx;
      for (int interval_idx = 0; interval_idx < channel_intervals->n_intervals;
           interval_idx++) {
        struct Interval const *const interval =
            channel_intervals->intervals + interval_idx;
        int n_times = interval->stop - interval->start;
        for (int time_idx = 0; time_idx < n_times; time_idx++) {
          sum_y += y[channel_idx][interval_idx][time_idx] *
                   yp[channel_idx][interval_idx][time_idx];
        }
      }
    }

    EXPECT_NEAR(sum_x, sum_y, 1e-12);

    free_model_fk(model.n_windows, &x);
    free_model_fk(model.n_windows, &xp);
    free_blended(&(model.blend_config), &y);
    free_blended(&(model.blend_config), &yp);
    free_model(&model);
  }
}
