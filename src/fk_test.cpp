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

TEST(FK, SetTaperWindow) {
  int window_length;
  int use_taper[2];
  AGD_TYPE *taper_window = nullptr;

  window_length = 2;
  use_taper[0] = 0;
  use_taper[1] = 0;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  free(taper_window);

  window_length = 2;
  use_taper[0] = 1;
  use_taper[1] = 1;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ZERO);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  free(taper_window);

  window_length = 2;
  use_taper[0] = 0;
  use_taper[1] = 1;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  free(taper_window);

  window_length = 2;
  use_taper[0] = 1;
  use_taper[1] = 0;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ZERO);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  free(taper_window);

  window_length = 3;
  use_taper[0] = 0;
  use_taper[1] = 0;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  free(taper_window);

  window_length = 3;
  use_taper[0] = 1;
  use_taper[1] = 1;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ZERO);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  free(taper_window);

  window_length = 3;
  use_taper[0] = 0;
  use_taper[1] = 1;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  free(taper_window);

  window_length = 3;
  use_taper[0] = 1;
  use_taper[1] = 0;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ZERO);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  free(taper_window);

  window_length = 4;
  use_taper[0] = 0;
  use_taper[1] = 0;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[3], AGD_ONE);
  free(taper_window);

  window_length = 4;
  use_taper[0] = 1;
  use_taper[1] = 1;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ZERO);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_HALF);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[3], AGD_HALF);
  free(taper_window);

  window_length = 4;
  use_taper[0] = 0;
  use_taper[1] = 1;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[3], AGD_HALF);
  free(taper_window);

  window_length = 4;
  use_taper[0] = 1;
  use_taper[1] = 0;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ZERO);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_HALF);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[3], AGD_ONE);
  free(taper_window);

  window_length = 5;
  use_taper[0] = 0;
  use_taper[1] = 0;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[3], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[4], AGD_ONE);
  free(taper_window);

  window_length = 5;
  use_taper[0] = 1;
  use_taper[1] = 1;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ZERO);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_HALF);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[3], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[4], AGD_HALF);
  free(taper_window);

  window_length = 5;
  use_taper[0] = 0;
  use_taper[1] = 1;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[3], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[4], AGD_HALF);
  free(taper_window);

  window_length = 5;
  use_taper[0] = 1;
  use_taper[1] = 0;
  EXPECT_EQ(set_taper_window(window_length, use_taper, &taper_window), 0);
  EXPECT_FLOAT_EQ(taper_window[0], AGD_ZERO);
  EXPECT_FLOAT_EQ(taper_window[1], AGD_HALF);
  EXPECT_FLOAT_EQ(taper_window[2], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[3], AGD_ONE);
  EXPECT_FLOAT_EQ(taper_window[4], AGD_ONE);
  free(taper_window);
}

TEST(FK, ForwardOnes) {
  int n_windows = 3;
  std::vector<int> n_times{16, 12, 16};
  std::vector<int> n_times_out_capacity{31, 16, 16};
  std::vector<int> n_dims{2, 2, 3};
  std::vector<std::vector<int>> window_shapes{{6, 8}, {5, 6}, {2, 3, 4}};
  std::vector<std::vector<int>> use_taper{{1, 0}, {1, 1}, {1, 1, 0, 0}};
  std::vector<int> window_fk_config_idx(n_windows);
  struct FKConfigs fk_configs = {};
  void **temp_1;
  void **temp_2;

  for (int window_idx = 0; window_idx < n_windows; window_idx++) {
    int n_times_window = window_shapes[window_idx][n_dims[window_idx] - 1];
    EXPECT_EQ(
        set_fk_configs(n_times[window_idx], n_times_out_capacity[window_idx],
                       n_times_window, n_dims[window_idx],
                       window_shapes[window_idx].data(),
                       use_taper[window_idx].data(), &fk_configs,
                       &window_fk_config_idx[window_idx]),
        0);
  }
  EXPECT_EQ(fk_configs.n_window_configs, 3);
  EXPECT_EQ(window_fk_config_idx[0], 0);
  EXPECT_EQ(window_fk_config_idx[1], 1);
  EXPECT_EQ(window_fk_config_idx[2], 2);
  EXPECT_EQ(fk_configs.window_configs[0].n_time_windows, 3);
  EXPECT_EQ(fk_configs.window_configs[0].n_dims, 2);
  EXPECT_EQ(fk_configs.window_configs[0].n_times, 16);
  EXPECT_EQ(fk_configs.window_configs[0].n_times_output_capacity, 31);
  EXPECT_EQ(fk_configs.window_configs[0].n_traces, 6);
  EXPECT_EQ(fk_configs.window_configs[0].n_traces_padded,
            (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[0].n_elements_tx, 6 * 16);
  EXPECT_EQ(fk_configs.window_configs[0].n_elements_tx_output_capacity, 6 * 31);
  EXPECT_EQ(fk_configs.window_configs[0].n_elements_padded_window,
            (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)) *
                (int)AGD_CEIL(8 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[0].n_elements_padded_window_fk,
            (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)) *
                ((int)AGD_CEIL(8 * (AGD_ONE + AGD_PAD_FRACTION)) / 2 + 1));
  EXPECT_EQ(fk_configs.window_configs[0].window_shape[0], 6);
  EXPECT_EQ(fk_configs.window_configs[0].window_shape[1], 8);
  EXPECT_EQ(fk_configs.window_configs[0].padded_window_shape[0],
            (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[0].padded_window_shape[1],
            (int)AGD_CEIL(8 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[0].padded_window_shape_fk[0],
            (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[0].padded_window_shape_fk[1],
            (int)AGD_CEIL(8 * (AGD_ONE + AGD_PAD_FRACTION)) / 2 + 1);
  EXPECT_EQ(fk_configs.window_configs[0].use_taper[0], 1);
  EXPECT_EQ(fk_configs.window_configs[0].use_taper[1], 0);
  EXPECT_FLOAT_EQ(
      fk_configs.window_configs[0].scaler,
      (AGD_TYPE)(AGD_ONE /
                 AGD_SQRT((int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)) *
                          (int)AGD_CEIL(8 * (AGD_ONE + AGD_PAD_FRACTION)))));
  for (int idx = 0; idx < 3; idx++) {
    EXPECT_FLOAT_EQ(
        fk_configs.window_configs[0].taper_windows[0][idx],
        AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI * (AGD_TYPE)idx /
                                      (AGD_TYPE)6)));
  }
  for (int idx = 3; idx < 6; idx++) {
    EXPECT_FLOAT_EQ(fk_configs.window_configs[0].taper_windows[0][idx],
                    AGD_ONE);
  }
  for (int idx = 0; idx < 8; idx++) {
    EXPECT_FLOAT_EQ(
        fk_configs.window_configs[0].taper_windows[1][idx],
        AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI * (AGD_TYPE)idx /
                                      (AGD_TYPE)8)));
  }

  EXPECT_EQ(fk_configs.window_configs[1].n_time_windows, 3);
  EXPECT_EQ(fk_configs.window_configs[1].n_dims, 2);
  EXPECT_EQ(fk_configs.window_configs[1].n_times, 12);
  EXPECT_EQ(fk_configs.window_configs[1].n_times_output_capacity, 16);
  EXPECT_EQ(fk_configs.window_configs[1].n_traces, 5);
  EXPECT_EQ(fk_configs.window_configs[1].n_traces_padded,
            (int)AGD_CEIL(5 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[1].n_elements_tx, 5 * 12);
  EXPECT_EQ(fk_configs.window_configs[1].n_elements_tx_output_capacity, 5 * 16);
  EXPECT_EQ(fk_configs.window_configs[1].n_elements_padded_window,
            (int)AGD_CEIL(5 * (AGD_ONE + AGD_PAD_FRACTION)) *
                (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[1].n_elements_padded_window_fk,
            (int)AGD_CEIL(5 * (AGD_ONE + AGD_PAD_FRACTION)) *
                ((int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)) / 2 + 1));
  EXPECT_EQ(fk_configs.window_configs[1].window_shape[0], 5);
  EXPECT_EQ(fk_configs.window_configs[1].window_shape[1], 6);
  EXPECT_EQ(fk_configs.window_configs[1].padded_window_shape[0],
            (int)AGD_CEIL(5 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[1].padded_window_shape[1],
            (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[1].padded_window_shape_fk[0],
            (int)AGD_CEIL(5 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[1].padded_window_shape_fk[1],
            (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)) / 2 + 1);
  EXPECT_EQ(fk_configs.window_configs[1].use_taper[0], 1);
  EXPECT_EQ(fk_configs.window_configs[1].use_taper[1], 1);
  EXPECT_FLOAT_EQ(
      fk_configs.window_configs[1].scaler,
      (AGD_TYPE)(AGD_ONE /
                 AGD_SQRT((int)AGD_CEIL(5 * (AGD_ONE + AGD_PAD_FRACTION)) *
                          (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)))));
  for (int idx = 0; idx < 2; idx++) {
    EXPECT_FLOAT_EQ(
        fk_configs.window_configs[1].taper_windows[0][idx],
        AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI * (AGD_TYPE)idx /
                                      (AGD_TYPE)4)));
  }
  EXPECT_FLOAT_EQ(fk_configs.window_configs[1].taper_windows[0][2], AGD_ONE);
  for (int idx = 0; idx < 2; idx++) {
    EXPECT_FLOAT_EQ(
        fk_configs.window_configs[1].taper_windows[0][3 + idx],
        AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                      (AGD_TYPE)(2 + idx) / (AGD_TYPE)4)));
  }
  for (int idx = 0; idx < 6; idx++) {
    EXPECT_FLOAT_EQ(
        fk_configs.window_configs[1].taper_windows[1][idx],
        AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI * (AGD_TYPE)idx /
                                      (AGD_TYPE)6)));
  }

  EXPECT_EQ(fk_configs.window_configs[2].n_time_windows, 7);
  EXPECT_EQ(fk_configs.window_configs[2].n_dims, 3);
  EXPECT_EQ(fk_configs.window_configs[2].n_times, 16);
  EXPECT_EQ(fk_configs.window_configs[2].n_times_output_capacity, 16);
  EXPECT_EQ(fk_configs.window_configs[2].n_traces, 6);
  EXPECT_EQ(fk_configs.window_configs[2].n_traces_padded,
            (int)AGD_CEIL(2 * (AGD_ONE + AGD_PAD_FRACTION)) *
                (int)AGD_CEIL(3 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[2].n_elements_tx, 2 * 3 * 16);
  EXPECT_EQ(fk_configs.window_configs[2].n_elements_tx_output_capacity,
            2 * 3 * 16);
  EXPECT_EQ(fk_configs.window_configs[2].n_elements_padded_window,
            (int)AGD_CEIL(2 * (AGD_ONE + AGD_PAD_FRACTION)) *
                (int)AGD_CEIL(3 * (AGD_ONE + AGD_PAD_FRACTION)) *
                (int)AGD_CEIL(4 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[2].n_elements_padded_window_fk,
            (int)AGD_CEIL(2 * (AGD_ONE + AGD_PAD_FRACTION)) *
                (int)AGD_CEIL(3 * (AGD_ONE + AGD_PAD_FRACTION)) *
                ((int)AGD_CEIL(4 * (AGD_ONE + AGD_PAD_FRACTION)) / 2 + 1));
  EXPECT_EQ(fk_configs.window_configs[2].window_shape[0], 2);
  EXPECT_EQ(fk_configs.window_configs[2].window_shape[1], 3);
  EXPECT_EQ(fk_configs.window_configs[2].window_shape[2], 4);
  EXPECT_EQ(fk_configs.window_configs[2].padded_window_shape[0],
            (int)AGD_CEIL(2 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[2].padded_window_shape[1],
            (int)AGD_CEIL(3 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[2].padded_window_shape[2],
            (int)AGD_CEIL(4 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[2].padded_window_shape_fk[0],
            (int)AGD_CEIL(2 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[2].padded_window_shape_fk[1],
            (int)AGD_CEIL(3 * (AGD_ONE + AGD_PAD_FRACTION)));
  EXPECT_EQ(fk_configs.window_configs[2].padded_window_shape_fk[2],
            (int)AGD_CEIL(4 * (AGD_ONE + AGD_PAD_FRACTION)) / 2 + 1);
  EXPECT_EQ(fk_configs.window_configs[2].use_taper[0], 1);
  EXPECT_EQ(fk_configs.window_configs[2].use_taper[1], 1);
  EXPECT_EQ(fk_configs.window_configs[2].use_taper[2], 0);
  EXPECT_EQ(fk_configs.window_configs[2].use_taper[3], 0);
  EXPECT_EQ(
      fk_configs.window_configs[2].scaler,
      (AGD_TYPE)(AGD_ONE /
                 AGD_SQRT((int)AGD_CEIL(2 * (AGD_ONE + AGD_PAD_FRACTION)) *
                          (int)AGD_CEIL(3 * (AGD_ONE + AGD_PAD_FRACTION)) *
                          (int)AGD_CEIL(4 * (AGD_ONE + AGD_PAD_FRACTION)))));
  for (int idx = 0; idx < 2; idx++) {
    EXPECT_FLOAT_EQ(
        fk_configs.window_configs[2].taper_windows[0][idx],
        AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI * (AGD_TYPE)idx /
                                      (AGD_TYPE)2)));
  }
  for (int idx = 0; idx < 3; idx++) {
    EXPECT_FLOAT_EQ(fk_configs.window_configs[2].taper_windows[1][idx],
                    AGD_ONE);
  }
  for (int idx = 0; idx < 4; idx++) {
    EXPECT_FLOAT_EQ(
        fk_configs.window_configs[2].taper_windows[2][idx],
        AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI * (AGD_TYPE)idx /
                                      (AGD_TYPE)4)));
  }

  EXPECT_EQ(allocate_temporary_arrays(&fk_configs, &temp_1, &temp_2), 0);
  EXPECT_EQ(set_fk_fftw_plans(temp_1[0], temp_2[0], &fk_configs), 0);

  int max_model_size = 0;
  for (int config_idx = 0; config_idx < fk_configs.n_window_configs;
       config_idx++) {
    struct FKWindowConfig const *const fk_config =
        fk_configs.window_configs + config_idx;
    int const model_size =
        fk_config->n_time_windows * fk_config->n_elements_padded_window_fk;
    if (model_size > max_model_size) max_model_size = model_size;
  }

  std::vector<AGD_FFTW_COMPLEX> model_fk(max_model_size);

  for (int window_idx = 0; window_idx < n_windows; window_idx++) {
    struct FKWindowConfig const *const fk_config =
        fk_configs.window_configs + window_fk_config_idx[window_idx];

    memset(model_fk.data(), 0, max_model_size * sizeof(AGD_FFTW_COMPLEX));
    for (int time_window_idx = 0; time_window_idx < fk_config->n_time_windows;
         time_window_idx++) {
      model_fk[time_window_idx * fk_config->n_elements_padded_window_fk][0] =
          AGD_ONE;
    }

    fk_forward(model_fk.data(), fk_config, temp_2[0], temp_1[0]);

    if (window_idx == 0) {
      AGD_TYPE scaler =
          AGD_ONE / AGD_SQRT((int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)) *
                             (int)AGD_CEIL(8 * (AGD_ONE + AGD_PAD_FRACTION)));
      for (int space0_idx = 0; space0_idx < 6; space0_idx++) {
        /* First three traces should be tapered */
        AGD_TYPE trace_scaler = AGD_ONE;
        if (space0_idx < 3) {
          trace_scaler =
              AGD_HALF *
              (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                 (AGD_TYPE)space0_idx / (AGD_TYPE)6));
        }

        /* First four time samples should be tapered */
        for (int time_idx = 0; time_idx < 4; time_idx++) {
          AGD_TYPE time_scaler =
              AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                            (AGD_TYPE)time_idx / (AGD_TYPE)8));
          EXPECT_FLOAT_EQ(((AGD_TYPE *)temp_1[0])[space0_idx * 31 + time_idx],
                          scaler * trace_scaler * time_scaler);
        }

        /* No taper on time samples 4 to 12 */
        for (int time_idx = 4; time_idx < 12; time_idx++) {
          EXPECT_FLOAT_EQ(((AGD_TYPE *)temp_1[0])[space0_idx * 31 + time_idx],
                          scaler * trace_scaler);
        }

        /* Last four time samples (in data portion of trace) should be tapered
         */
        for (int time_idx = 12; time_idx < 16; time_idx++) {
          AGD_TYPE time_scaler =
              AGD_HALF *
              (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                 (AGD_TYPE)(time_idx + 4 - 12) / (AGD_TYPE)8));
          EXPECT_FLOAT_EQ(((AGD_TYPE *)temp_1[0])[space0_idx * 31 + time_idx],
                          scaler * trace_scaler * time_scaler);
        }

        /* The trace should be zero from the end of the data portion */
        for (int time_idx = 16; time_idx < 31; time_idx++) {
          EXPECT_FLOAT_EQ(((AGD_TYPE *)temp_1[0])[space0_idx * 31 + time_idx],
                          AGD_ZERO);
        }
      }
    } else if (window_idx == 1) {
      AGD_TYPE scaler =
          AGD_ONE / AGD_SQRT((int)AGD_CEIL(5 * (AGD_ONE + AGD_PAD_FRACTION)) *
                             (int)AGD_CEIL(6 * (AGD_ONE + AGD_PAD_FRACTION)));
      for (int space0_idx = 0; space0_idx < 5; space0_idx++) {
        /* First two and last two traces should be tapered */
        AGD_TYPE trace_scaler = AGD_ONE;
        if (space0_idx < 2) {
          trace_scaler =
              AGD_HALF *
              (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                 (AGD_TYPE)space0_idx / (AGD_TYPE)4));
        }
        if (space0_idx > 2) {
          trace_scaler =
              AGD_HALF *
              (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                 (AGD_TYPE)(space0_idx - 1) / (AGD_TYPE)4));
        }

        /* First three time samples should be tapered */
        for (int time_idx = 0; time_idx < 3; time_idx++) {
          AGD_TYPE time_scaler =
              AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                            (AGD_TYPE)time_idx / (AGD_TYPE)6));
          EXPECT_FLOAT_EQ(((AGD_TYPE *)temp_1[0])[space0_idx * 16 + time_idx],
                          scaler * trace_scaler * time_scaler);
        }

        /* No taper on time samples 3 to 9 */
        for (int time_idx = 3; time_idx < 9; time_idx++) {
          EXPECT_FLOAT_EQ(((AGD_TYPE *)temp_1[0])[space0_idx * 16 + time_idx],
                          scaler * trace_scaler);
        }

        /* Last three time samples (in data portion of trace) should be tapered
         */
        for (int time_idx = 9; time_idx < 12; time_idx++) {
          AGD_TYPE time_scaler =
              AGD_HALF *
              (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                 (AGD_TYPE)(time_idx + 3 - 9) / (AGD_TYPE)6));
          EXPECT_FLOAT_EQ(((AGD_TYPE *)temp_1[0])[space0_idx * 16 + time_idx],
                          scaler * trace_scaler * time_scaler);
        }

        /* The trace should be zero from the end of the data portion */
        for (int time_idx = 12; time_idx < 16; time_idx++) {
          EXPECT_FLOAT_EQ(((AGD_TYPE *)temp_1[0])[space0_idx * 16 + time_idx],
                          AGD_ZERO);
        }
      }
    } else {
      AGD_TYPE scaler =
          AGD_ONE / AGD_SQRT((int)AGD_CEIL(2 * (AGD_ONE + AGD_PAD_FRACTION)) *
                             (int)AGD_CEIL(3 * (AGD_ONE + AGD_PAD_FRACTION)) *
                             (int)AGD_CEIL(4 * (AGD_ONE + AGD_PAD_FRACTION)));

      /* The first trace of the first dimension should be zero */
      for (int space1_idx = 0; space1_idx < 3; space1_idx++) {
        for (int time_idx = 0; time_idx < 16; time_idx++) {
          EXPECT_FLOAT_EQ(((AGD_TYPE *)temp_1[0])[space1_idx * 16 + time_idx],
                          AGD_ZERO);
        }
      }

      /* The second (final) trace of the first dimension: */

      /* First two time samples should be tapered */
      for (int space1_idx = 0; space1_idx < 3; space1_idx++) {
        for (int time_idx = 0; time_idx < 2; time_idx++) {
          AGD_TYPE time_scaler =
              AGD_HALF * (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                            (AGD_TYPE)time_idx / (AGD_TYPE)4));
          EXPECT_FLOAT_EQ(
              ((AGD_TYPE *)temp_1[0])[1 * 3 * 16 + space1_idx * 16 + time_idx],
              scaler * time_scaler);
        }
      }

      /* Time samples 1 to 15 should have no taper */
      for (int space1_idx = 0; space1_idx < 3; space1_idx++) {
        for (int time_idx = 2; time_idx < 14; time_idx++) {
          EXPECT_FLOAT_EQ(
              ((AGD_TYPE *)temp_1[0])[1 * 3 * 16 + space1_idx * 16 + time_idx],
              scaler);
        }
      }

      /* Last two time samples should be tapered */
      for (int space1_idx = 0; space1_idx < 3; space1_idx++) {
        for (int time_idx = 14; time_idx < 16; time_idx++) {
          AGD_TYPE time_scaler =
              AGD_HALF *
              (AGD_ONE - AGD_COS(AGD_TWO * (AGD_TYPE)M_PI *
                                 (AGD_TYPE)(time_idx + 2 - 14) / (AGD_TYPE)4));
          EXPECT_FLOAT_EQ(
              ((AGD_TYPE *)temp_1[0])[1 * 3 * 16 + space1_idx * 16 + time_idx],
              scaler * time_scaler);
        }
      }
    }
  }

  free_temporary_arrays(&temp_1, &temp_2);
  free_fk_config(&fk_configs);
}

TEST(FK, DotTest) {
  srand(1);
  for (int rand_repeat = 0; rand_repeat < 5; rand_repeat++) {
    int n_windows = 1 + rand() % 4;
    std::vector<int> n_times(n_windows);
    std::vector<int> n_times_out_capacity(n_windows);
    std::vector<int> n_times_window(n_windows);
    std::vector<int> n_dims(n_windows);
    std::vector<std::vector<int>> window_shapes(n_windows);
    std::vector<std::vector<int>> use_taper(n_windows);
    std::vector<int> window_fk_config_idx(n_windows);
    struct FKConfigs fk_configs = {};
    void **temp_1;
    void **temp_2;

    for (int window_idx = 0; window_idx < n_windows; window_idx++) {
      n_dims[window_idx] = 2 + rand() % 4;
      for (int dim_idx = 0; dim_idx < n_dims[window_idx]; dim_idx++) {
        if (dim_idx == n_dims[window_idx] - 1) {
          window_shapes[window_idx].push_back(2 + 2 * (rand() % 8));
        } else {
          window_shapes[window_idx].push_back(2 + rand() % 8);
        }
        /* no edges in time dimension */
        if (dim_idx == n_dims[window_idx] - 1) {
          use_taper[window_idx].push_back(1);
          use_taper[window_idx].push_back(1);
        } else {
          use_taper[window_idx].push_back(rand() % 2);
          use_taper[window_idx].push_back(rand() % 2);
        }
      }
      n_times_window[window_idx] =
          window_shapes[window_idx][n_dims[window_idx] - 1];
      n_times[window_idx] = n_times_window[window_idx] / 2 * (2 + rand() % 3);
      n_times_out_capacity[window_idx] = n_times[window_idx] + rand() % 16;

      EXPECT_EQ(
          set_fk_configs(n_times[window_idx], n_times_out_capacity[window_idx],
                         n_times_window[window_idx], n_dims[window_idx],
                         window_shapes[window_idx].data(),
                         use_taper[window_idx].data(), &fk_configs,
                         &window_fk_config_idx[window_idx]),
          0);
    }

    EXPECT_EQ(allocate_temporary_arrays(&fk_configs, &temp_1, &temp_2), 0);
    EXPECT_EQ(set_fk_fftw_plans(temp_1[0], temp_2[0], &fk_configs), 0);

    for (int window_idx = 0; window_idx < n_windows; window_idx++) {
      struct FKWindowConfig const *const fk_config =
          fk_configs.window_configs + window_fk_config_idx[window_idx];
      int const model_size =
          fk_config->n_time_windows * fk_config->n_elements_padded_window_fk;
      int const out_size = fk_config->n_elements_tx_output_capacity;

      AGD_FFTW_COMPLEX *x;
      AGD_FFTW_COMPLEX *xp;
      AGD_TYPE *y;
      AGD_TYPE *yp;

      x = (AGD_FFTW_COMPLEX *)AGD_FFTW_ALLOC_COMPLEX(model_size);
      xp = (AGD_FFTW_COMPLEX *)AGD_FFTW_ALLOC_COMPLEX(model_size);
      y = (AGD_TYPE *)AGD_FFTW_MALLOC((size_t)out_size * sizeof(AGD_TYPE));
      yp = (AGD_TYPE *)AGD_FFTW_MALLOC((size_t)out_size * sizeof(AGD_TYPE));

      /* Make x by transforming random real numbers */
      for (int idx = 0; idx < fk_config->n_time_windows *
                                  fk_config->n_elements_padded_window;
           idx++) {
        static_cast<AGD_TYPE *>(temp_2[0])[idx] =
            (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX - AGD_HALF;
      }
      AGD_FFTW_EXECUTE(fk_config->adj_plan);
      memcpy(x, temp_1[0], model_size * sizeof(AGD_FFTW_COMPLEX));

      /* Make yp */
      for (int idx = 0; idx < out_size; idx++) {
        yp[idx] = (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX - AGD_HALF;
      }

      /* x -> y */
      fk_forward(x, fk_config, temp_2[0], temp_1[0]);
      memcpy(y, temp_1[0], out_size * sizeof(AGD_TYPE));

      /* yp -> xp */
      memcpy(temp_1[0], yp, out_size * sizeof(AGD_TYPE));
      fk_adjoint(fk_config, temp_2[0], temp_1[0]);
      memcpy(xp, temp_1[0], model_size * sizeof(AGD_FFTW_COMPLEX));

      /* x sum */
      AGD_TYPE sum_x = AGD_ZERO;
      for (int idx = 0; idx < model_size; idx++) {
        sum_x += x[idx][0] * xp[idx][0] + x[idx][1] * xp[idx][1];
      }

      /* y sum */
      AGD_TYPE sum_y = AGD_ZERO;
      for (int idx = 0; idx < out_size; idx++) {
        sum_y += y[idx] * yp[idx];
      }

      EXPECT_NEAR(sum_x, sum_y, 5e-11);

      AGD_FFTW_FREE(x);
      AGD_FFTW_FREE(xp);
      AGD_FFTW_FREE(y);
      AGD_FFTW_FREE(yp);
    }

    free_temporary_arrays(&temp_1, &temp_2);
    free_fk_config(&fk_configs);
  }
}
