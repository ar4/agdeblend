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

TEST(Wavelet, Forward) {
  int n_traces[2] = {2, 4};
  int n_times[2] = {3 + 2, 4 + 2};
  std::vector<std::vector<enum AGDTraceType>> trace_types = {
      {AGDLive, AGDLive}, {AGDLive, AGDLive, AGDLive, AGDMissing}};
  std::vector<std::vector<int>> wavelet_idxs = {{0, 1}, {0, 1, 2, 50}};
  struct WaveletConfig config = {};
  std::vector<struct WaveletParams> params(2);
  std::vector<std::vector<AGD_TYPE>> wavelets = {
      {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  std::vector<AGD_TYPE *> wavelets_ptr = {&wavelets[0][0], &wavelets[1][0],
                                          &wavelets[2][0]};
  std::vector<int> wavelet_lengths = {3, 3, 3};
  std::vector<std::vector<AGD_TYPE>> data(2);

  data[0].resize(n_traces[0] * n_times[0]);
  std::iota(data[0].begin(), data[0].end(), 5);
  data[1].resize(n_traces[1] * n_times[1]);
  std::iota(data[1].begin(), data[1].end(), 100);

  for (int idx = 0; idx < 2; idx++) {
    EXPECT_EQ(
        set_wavelet_params(n_traces[idx], n_times[idx], trace_types[idx].data(),
                           wavelet_idxs[idx].data(), &config, &params[idx]),
        0);
  }

  EXPECT_EQ(config.n_source_configs, 5);
  EXPECT_EQ(config.n_fftw_plans, 2);
  EXPECT_EQ(config.source_configs[0].original_idx, 0);
  EXPECT_EQ(config.source_configs[0].n_times, n_times[0]);
  EXPECT_EQ(config.source_configs[1].original_idx, 1);
  EXPECT_EQ(config.source_configs[1].n_times, n_times[0]);
  EXPECT_EQ(config.source_configs[2].original_idx, 0);
  EXPECT_EQ(config.source_configs[2].n_times, n_times[1]);
  EXPECT_EQ(config.source_configs[3].original_idx, 1);
  EXPECT_EQ(config.source_configs[3].n_times, n_times[1]);
  EXPECT_EQ(config.source_configs[4].original_idx, 2);
  EXPECT_EQ(config.source_configs[4].n_times, n_times[1]);
  EXPECT_EQ(config.plans[0].n_traces, 2);
  EXPECT_EQ(config.plans[0].n_times, n_times[0]);
  EXPECT_EQ(config.plans[1].n_traces, 4);
  EXPECT_EQ(config.plans[1].n_times, n_times[1]);
  EXPECT_EQ(params[0].plan_idx, 0);
  EXPECT_EQ(params[0].source_config_idxs[0], 0);
  EXPECT_EQ(params[0].source_config_idxs[1], 1);
  EXPECT_EQ(params[1].plan_idx, 1);
  EXPECT_EQ(params[1].source_config_idxs[0], 2);
  EXPECT_EQ(params[1].source_config_idxs[1], 3);
  EXPECT_EQ(params[1].source_config_idxs[2], 4);
  EXPECT_EQ(params[1].source_config_idxs[3], -1);

  std::vector<AGD_TYPE> temp_1(2 * n_traces[1] * (n_times[1] / 2 + 1));
  std::vector<AGD_TYPE> temp_2(2 * n_traces[1] * (n_times[1] / 2 + 1));

  EXPECT_EQ(set_wavelet_fftw_plans(temp_1.data(), temp_2.data(), &config), 0);
  EXPECT_EQ(transform_wavelets(&config, wavelet_lengths.data(),
                               wavelets_ptr.data(), temp_1.data()),
            0);

  const AGD_TYPE eps = 1e-12;

  for (int trace_idx = 0; trace_idx < n_traces[0]; trace_idx++) {
    for (int time_idx = 0; time_idx < 3; time_idx++) {
      temp_1[trace_idx * 5 + time_idx] = data[0][trace_idx * 5 + time_idx];
    }
  }
  wavelet_forward(&config, &params[0], temp_1.data(), temp_2.data());
  EXPECT_NEAR(temp_1[0], (AGD_TYPE)5, eps);
  EXPECT_NEAR(temp_1[1], (AGD_TYPE)6, eps);
  EXPECT_NEAR(temp_1[2], (AGD_TYPE)7, eps);
  EXPECT_NEAR(temp_1[3], (AGD_TYPE)0, eps);
  EXPECT_NEAR(temp_1[4], (AGD_TYPE)0, eps);

  EXPECT_NEAR(temp_1[5], (AGD_TYPE)0, eps);
  EXPECT_NEAR(temp_1[6], (AGD_TYPE)10, eps);
  EXPECT_NEAR(temp_1[7], (AGD_TYPE)11, eps);
  EXPECT_NEAR(temp_1[8], (AGD_TYPE)12, eps);
  EXPECT_NEAR(temp_1[9], (AGD_TYPE)0, eps);

  for (int trace_idx = 0; trace_idx < n_traces[1]; trace_idx++) {
    for (int time_idx = 0; time_idx < 4; time_idx++) {
      temp_1[trace_idx * 6 + time_idx] = data[1][trace_idx * 6 + time_idx];
    }
  }
  wavelet_forward(&config, &params[1], temp_1.data(), temp_2.data());
  EXPECT_NEAR(temp_1[0], (AGD_TYPE)100, eps);
  EXPECT_NEAR(temp_1[1], (AGD_TYPE)101, eps);
  EXPECT_NEAR(temp_1[2], (AGD_TYPE)102, eps);
  EXPECT_NEAR(temp_1[3], (AGD_TYPE)103, eps);
  EXPECT_NEAR(temp_1[4], (AGD_TYPE)0, eps);
  EXPECT_NEAR(temp_1[5], (AGD_TYPE)0, eps);

  EXPECT_NEAR(temp_1[6], (AGD_TYPE)0, eps);
  EXPECT_NEAR(temp_1[7], (AGD_TYPE)106, eps);
  EXPECT_NEAR(temp_1[8], (AGD_TYPE)107, eps);
  EXPECT_NEAR(temp_1[9], (AGD_TYPE)108, eps);
  EXPECT_NEAR(temp_1[10], (AGD_TYPE)109, eps);
  EXPECT_NEAR(temp_1[11], (AGD_TYPE)0, eps);

  EXPECT_NEAR(temp_1[12], (AGD_TYPE)0, eps);
  EXPECT_NEAR(temp_1[13], (AGD_TYPE)0, eps);
  EXPECT_NEAR(temp_1[14], (AGD_TYPE)112, eps);
  EXPECT_NEAR(temp_1[15], (AGD_TYPE)113, eps);
  EXPECT_NEAR(temp_1[16], (AGD_TYPE)114, eps);
  EXPECT_NEAR(temp_1[17], (AGD_TYPE)115, eps);

  for (int idx = 0; idx < 2; idx++) {
    free_wavelet_params(&params[idx]);
  }
  free_wavelet_config(&config);
}

TEST(Wavelet, DotTest) {
  srand(1);
  for (int rand_repeat = 0; rand_repeat < 5; rand_repeat++) {
    int n_volumes = rand() % 5;
    std::vector<int> n_traces(n_volumes);
    std::vector<int> n_times_preconv(n_volumes);
    std::vector<int> n_times_conv(n_volumes);
    std::vector<std::vector<enum AGDTraceType>> trace_types(n_volumes);
    std::vector<enum AGDTraceType *> trace_types_ptr(n_volumes, nullptr);
    std::vector<std::vector<int>> wavelet_idxs(n_volumes);
    struct WaveletConfig config = {};
    std::vector<struct WaveletParams> params(n_volumes);
    std::vector<std::vector<AGD_TYPE>> wavelets;
    std::vector<AGD_TYPE *> wavelets_ptr;
    std::vector<int> wavelet_lengths;
    std::vector<std::vector<AGD_TYPE>> x(n_volumes);
    std::vector<std::vector<AGD_TYPE>> xp(n_volumes);
    std::vector<std::vector<AGD_TYPE>> y(n_volumes);
    std::vector<std::vector<AGD_TYPE>> yp(n_volumes);

    for (int volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      n_traces[volume_idx] = 2 + rand() % 8;
      n_times_preconv[volume_idx] = 1 + rand() % 8;
      n_times_conv[volume_idx] = n_times_preconv[volume_idx] + rand() % 4;
      trace_types[volume_idx].resize(n_traces[volume_idx]);
      std::generate(trace_types[volume_idx].begin(),
                    trace_types[volume_idx].end(),
                    []() { return (enum AGDTraceType)(rand() % 3); });

      int wavelet_length =
          n_times_conv[volume_idx] - n_times_preconv[volume_idx] + 1;
      /* Create some new wavelets of this length */
      int n_wavelets_to_create = 1 + rand() % (n_traces[volume_idx] - 1);
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
      for (int trace_idx = 0; trace_idx < n_traces[volume_idx]; trace_idx++) {
        wavelet_idxs[volume_idx].push_back(possible_idxs[rand() % n_possible]);
      }

      /* Create x, xp */
      x[volume_idx].resize(n_traces[volume_idx] * n_times_preconv[volume_idx]);
      xp[volume_idx].resize(n_traces[volume_idx] * n_times_preconv[volume_idx]);
      std::generate(x[volume_idx].begin(), x[volume_idx].end(), [] {
        return (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX - AGD_HALF;
      });

      /* Create y, yp */
      y[volume_idx].resize(n_traces[volume_idx] * n_times_conv[volume_idx]);
      yp[volume_idx].resize(n_traces[volume_idx] * n_times_conv[volume_idx]);
      std::generate(yp[volume_idx].begin(), yp[volume_idx].end(), [] {
        return (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX - AGD_HALF;
      });
    }

    /* Set wavelets_ptr */
    for (int wavelet_idx = 0;
         static_cast<size_t>(wavelet_idx) < wavelets.size(); wavelet_idx++) {
      wavelets_ptr.push_back(&wavelets[wavelet_idx][0]);
    }

    /* Params */
    for (int volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      EXPECT_EQ(
          set_wavelet_params(n_traces[volume_idx], n_times_conv[volume_idx],
                             trace_types[volume_idx].data(),
                             wavelet_idxs[volume_idx].data(), &config,
                             &params[volume_idx]),
          0);
    }

    /* FFTW plans and wavelet transformation */
    int n_workspace = 0;
    for (int volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      int volume_workspace =
          2 * n_traces[volume_idx] * (n_times_conv[volume_idx] / 2 + 1);
      if (volume_workspace > n_workspace) n_workspace = volume_workspace;
    }
    std::vector<AGD_TYPE> temp_1(n_workspace);
    std::vector<AGD_TYPE> temp_2(n_workspace);
    EXPECT_EQ(set_wavelet_fftw_plans(temp_1.data(), temp_2.data(), &config), 0);
    EXPECT_EQ(transform_wavelets(&config, wavelet_lengths.data(),
                                 wavelets_ptr.data(), temp_1.data()),
              0);

    /* Forward */
    for (int volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      memset(
          temp_1.data(), 0,
          n_traces[volume_idx] * n_times_conv[volume_idx] * sizeof(AGD_TYPE));
      memcpy(temp_1.data(), x[volume_idx].data(),
             n_traces[volume_idx] * n_times_preconv[volume_idx] *
                 sizeof(AGD_TYPE));
      wavelet_forward(&config, &params[volume_idx], temp_1.data(),
                      temp_2.data());
      memcpy(
          y[volume_idx].data(), temp_1.data(),
          n_traces[volume_idx] * n_times_conv[volume_idx] * sizeof(AGD_TYPE));
    }

    /* Adjoint */
    for (int volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      memset(
          temp_1.data(), 0,
          n_traces[volume_idx] * n_times_conv[volume_idx] * sizeof(AGD_TYPE));
      memcpy(
          temp_1.data(), yp[volume_idx].data(),
          n_traces[volume_idx] * n_times_conv[volume_idx] * sizeof(AGD_TYPE));
      wavelet_adjoint(&config, &params[volume_idx], temp_1.data(),
                      temp_2.data());
      memcpy(xp[volume_idx].data(), temp_1.data(),
             n_traces[volume_idx] * n_times_preconv[volume_idx] *
                 sizeof(AGD_TYPE));
    }

    /* Sum x */
    AGD_TYPE sum_x = AGD_ZERO;
    for (int volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      sum_x += std::inner_product(x[volume_idx].cbegin(), x[volume_idx].cend(),
                                  xp[volume_idx].begin(), AGD_ZERO);
    }

    /* Sum y */
    AGD_TYPE sum_y = AGD_ZERO;
    for (int volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      sum_y += std::inner_product(y[volume_idx].cbegin(), y[volume_idx].cend(),
                                  yp[volume_idx].begin(), AGD_ZERO);
    }

    EXPECT_NEAR(sum_x, sum_y, 1e-12);

    for (int volume_idx = 0; volume_idx < n_volumes; volume_idx++) {
      free_wavelet_params(&params[volume_idx]);
    }
    free_wavelet_config(&config);
  }
}
