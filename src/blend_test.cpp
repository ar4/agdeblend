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

TEST(AddInterval, FirstInterval) {
  struct ChannelIntervals channel_intervals = {};
  struct Interval interval;
  interval.start = 100;
  interval.stop = 200;

  EXPECT_EQ(add_interval(&interval, &channel_intervals), 0);
  EXPECT_EQ(channel_intervals.n_intervals, 1);
  EXPECT_EQ(channel_intervals.intervals[0].start, 100);
  EXPECT_EQ(channel_intervals.intervals[0].stop, 200);
  free(channel_intervals.intervals);
}

TEST(AddInterval, TwoDisjointIntervals) {
  struct ChannelIntervals channel_intervals = {};
  struct Interval interval1;
  struct Interval interval2;
  interval1.start = 100;
  interval1.stop = 200;
  interval2.start = 300;
  interval2.stop = 400;

  EXPECT_EQ(add_interval(&interval1, &channel_intervals), 0);
  EXPECT_EQ(add_interval(&interval2, &channel_intervals), 0);
  EXPECT_EQ(channel_intervals.n_intervals, 2);
  EXPECT_EQ(channel_intervals.intervals[0].start, 100);
  EXPECT_EQ(channel_intervals.intervals[0].stop, 200);
  EXPECT_EQ(channel_intervals.intervals[1].start, 300);
  EXPECT_EQ(channel_intervals.intervals[1].stop, 400);
  free(channel_intervals.intervals);
}

TEST(AddInterval, TwoOverlappingIntervals) {
  struct ChannelIntervals channel_intervals = {};
  struct Interval interval1;
  struct Interval interval2;
  interval1.start = 100;
  interval1.stop = 200;
  interval2.start = 150;
  interval2.stop = 250;

  EXPECT_EQ(add_interval(&interval1, &channel_intervals), 0);
  EXPECT_EQ(add_interval(&interval2, &channel_intervals), 0);
  EXPECT_EQ(channel_intervals.n_intervals, 1);
  EXPECT_EQ(channel_intervals.intervals[0].start, 100);
  EXPECT_EQ(channel_intervals.intervals[0].stop, 250);
  free(channel_intervals.intervals);
}

TEST(AddInterval, ThreeOverlappingIntervals) {
  struct ChannelIntervals channel_intervals = {};
  struct Interval interval1;
  struct Interval interval2;
  struct Interval interval3;
  struct Interval interval4;
  struct Interval interval5;
  interval1.start = 50;
  interval1.stop = 75;
  interval2.start = 100;
  interval2.stop = 200;
  interval3.start = 250;
  interval3.stop = 350;
  interval4.start = 400;
  interval4.stop = 500;
  interval5.start = 175;
  interval5.stop = 275;

  EXPECT_EQ(add_interval(&interval1, &channel_intervals), 0);
  EXPECT_EQ(add_interval(&interval2, &channel_intervals), 0);
  EXPECT_EQ(add_interval(&interval3, &channel_intervals), 0);
  EXPECT_EQ(add_interval(&interval4, &channel_intervals), 0);
  EXPECT_EQ(add_interval(&interval5, &channel_intervals), 0);
  EXPECT_EQ(channel_intervals.n_intervals, 3);
  EXPECT_EQ(channel_intervals.intervals[0].start, 50);
  EXPECT_EQ(channel_intervals.intervals[0].stop, 75);
  EXPECT_EQ(channel_intervals.intervals[1].start, 100);
  EXPECT_EQ(channel_intervals.intervals[1].stop, 350);
  EXPECT_EQ(channel_intervals.intervals[2].start, 400);
  EXPECT_EQ(channel_intervals.intervals[2].stop, 500);
  free(channel_intervals.intervals);
}

TEST(GetUniqueChannels, TwoPatchs) {
  int n_patches = 2;
  std::vector<int> n_channels{4, 2};
  std::vector<std::vector<int>> channels{std::vector<int>(n_channels[0]),
                                         std::vector<int>(n_channels[1])};
  std::vector<int *> channels_ptr(n_patches, nullptr);
  std::vector<std::vector<enum AGDTraceType>> trace_types{
      std::vector<AGDTraceType>(n_channels[0], AGDLive),
      std::vector<AGDTraceType>(n_channels[1], AGDLive)};
  std::vector<AGDTraceType *> trace_types_ptr(n_patches, nullptr);
  int *unique_channel_values = nullptr;
  int n_unique_channels;

  channels[0][0] = 123;
  channels[0][1] = 456;
  channels[0][2] = 0;
  channels[0][3] = 123;
  channels[1][0] = 123;
  channels[1][1] = -789;

  trace_types[0][1] = AGDBad;
  trace_types[0][2] = AGDMissing;

  for (int i = 0; i < n_patches; i++) {
    channels_ptr[i] = &channels[i][0];
    trace_types_ptr[i] = &trace_types[i][0];
  }

  EXPECT_EQ(
      set_unique_channel_values(n_patches, n_channels.data(),
                                channels_ptr.data(), trace_types_ptr.data(),
                                &n_unique_channels, &unique_channel_values),
      0);

  EXPECT_EQ(n_unique_channels, 3);
  EXPECT_EQ(unique_channel_values[0], 123);
  EXPECT_EQ(unique_channel_values[1], 456);
  EXPECT_EQ(unique_channel_values[2], -789);
  EXPECT_EQ(get_channel_idx(channels[0][0], trace_types[0][0],
                            n_unique_channels, unique_channel_values),
            0);
  EXPECT_EQ(get_channel_idx(channels[0][1], trace_types[0][1],
                            n_unique_channels, unique_channel_values),
            1);
  EXPECT_EQ(get_channel_idx(channels[0][2], trace_types[0][2],
                            n_unique_channels, unique_channel_values),
            -1);
  EXPECT_EQ(get_channel_idx(channels[0][3], trace_types[0][3],
                            n_unique_channels, unique_channel_values),
            0);
  EXPECT_EQ(get_channel_idx(channels[1][0], trace_types[1][0],
                            n_unique_channels, unique_channel_values),
            0);
  EXPECT_EQ(get_channel_idx(channels[1][1], trace_types[1][1],
                            n_unique_channels, unique_channel_values),
            2);

  free(unique_channel_values);
}

TEST(GetBlendConfig, TwoChannels) {
  int n_patches = 2;
  int n_traces_per_patch[] = {1, 1};
  int n_times_per_patch[] = {100, 50};
  long int **shottimes;
  int **channels;
  enum AGDTraceType **trace_types;
  struct BlendConfig blend_config = {};

  shottimes = (long int **)malloc(2 * sizeof(long int *));
  channels = (int **)malloc(2 * sizeof(int *));
  trace_types = (AGDTraceType **)malloc(2 * sizeof(AGDTraceType *));
  shottimes[0] = (long int *)malloc(sizeof(long int));
  shottimes[1] = (long int *)malloc(sizeof(long int));
  channels[0] = (int *)malloc(sizeof(int));
  channels[1] = (int *)malloc(sizeof(int));
  trace_types[0] = (AGDTraceType *)malloc(sizeof(AGDTraceType));
  trace_types[1] = (AGDTraceType *)malloc(sizeof(AGDTraceType));

  shottimes[0][0] = 25;
  shottimes[1][0] = 35;
  channels[0][0] = 1;
  channels[1][0] = 0;
  trace_types[0][0] = AGDLive;
  trace_types[1][0] = AGDLive;

  EXPECT_EQ(set_blend_config(n_patches, n_traces_per_patch, n_times_per_patch,
                             shottimes, channels, trace_types, &blend_config),
            0);
  EXPECT_EQ(blend_config.n_channels, 2);
  EXPECT_EQ(blend_config.unique_channel_values[0], 1);
  EXPECT_EQ(blend_config.unique_channel_values[1], 0);
  EXPECT_EQ(blend_config.channels_intervals[0].n_intervals, 1);
  EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].start, 25);
  EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].stop, 25 + 100);
  EXPECT_EQ(blend_config.channels_intervals[1].n_intervals, 1);
  EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].start, 35);
  EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].stop, 35 + 50);

  free(shottimes[0]);
  free(shottimes[1]);
  free(shottimes);
  free(channels[0]);
  free(channels[1]);
  free(channels);
  free(trace_types[0]);
  free(trace_types[1]);
  free(trace_types);
  free_blend_config(&blend_config);
}

TEST(GetBlendParams, TwoLive) {
  int n_traces = 2;
  int n_times = 50;
  long int shottimes[] = {100, 200};
  int channels[] = {0, 0};
  enum AGDTraceType trace_types[] = {AGDLive, AGDLive};
  struct Interval *intervals;
  struct ChannelIntervals *channels_intervals;
  struct BlendConfig blend_config = {};
  struct BlendParams blend_params = {};

  intervals = (Interval *)malloc(sizeof(Interval));
  channels_intervals = (ChannelIntervals *)malloc(sizeof(ChannelIntervals));
  intervals[0].start = 100;
  intervals[0].stop = 250;
  channels_intervals[0].n_intervals = 1;
  channels_intervals[0].intervals = intervals;
  blend_config.n_channels = 1;
  blend_config.unique_channel_values = (int *)malloc(sizeof(int));
  blend_config.unique_channel_values[0] = 0;
  blend_config.channels_intervals = channels_intervals;

  EXPECT_EQ(set_blend_params(n_traces, n_times, shottimes, channels,
                             trace_types, &blend_config, &blend_params),
            0);
  EXPECT_EQ(blend_params.n_traces, 2);
  EXPECT_EQ(blend_params.trace_length, 50);
  EXPECT_EQ(blend_params.blend_coords[0].channel_idx, 0);
  EXPECT_EQ(blend_params.blend_coords[0].interval_idx, 0);
  EXPECT_EQ(blend_params.blend_coords[0].interval_start, 0);
  EXPECT_EQ(blend_params.blend_coords[0].trace_start, 0);
  EXPECT_EQ(blend_params.blend_coords[0].n_times, 50);
  EXPECT_EQ(blend_params.blend_coords[1].channel_idx, 0);
  EXPECT_EQ(blend_params.blend_coords[1].interval_idx, 0);
  EXPECT_EQ(blend_params.blend_coords[1].interval_start, 100);
  EXPECT_EQ(blend_params.blend_coords[1].trace_start, 0);
  EXPECT_EQ(blend_params.blend_coords[1].n_times, 50);

  free_blend_params(&blend_params);
  free_blend_config(&blend_config);
}

TEST(GetBlendParams, PartialOverlap) {
  int n_traces = 1;
  int n_times = 100;
  long int shottimes[] = {100};
  int channels[] = {0};
  enum AGDTraceType trace_types[] = {AGDLive};
  struct Interval *intervals;
  struct ChannelIntervals *channels_intervals;
  struct BlendConfig blend_config = {};
  struct BlendParams blend_params = {};

  intervals = (Interval *)malloc(sizeof(Interval));
  channels_intervals = (ChannelIntervals *)malloc(sizeof(ChannelIntervals));
  intervals[0].start = 50;
  intervals[0].stop = 150;
  channels_intervals[0].n_intervals = 1;
  channels_intervals[0].intervals = intervals;
  blend_config.n_channels = 1;
  blend_config.unique_channel_values = (int *)malloc(sizeof(int));
  blend_config.unique_channel_values[0] = 0;
  blend_config.channels_intervals = channels_intervals;

  EXPECT_EQ(set_blend_params(n_traces, n_times, shottimes, channels,
                             trace_types, &blend_config, &blend_params),
            0);
  EXPECT_EQ(blend_params.n_traces, 1);
  EXPECT_EQ(blend_params.trace_length, 100);
  EXPECT_EQ(blend_params.blend_coords[0].channel_idx, 0);
  EXPECT_EQ(blend_params.blend_coords[0].interval_idx, 0);
  EXPECT_EQ(blend_params.blend_coords[0].interval_start, 50);
  EXPECT_EQ(blend_params.blend_coords[0].trace_start, 0);
  EXPECT_EQ(blend_params.blend_coords[0].n_times, 50);

  free_blend_params(&blend_params);
  free_blend_config(&blend_config);
}

TEST(GetBlendParams, Missing) {
  int n_traces = 1;
  int n_times = 50;
  long int shottimes[] = {-1};
  int channels[] = {-1};
  enum AGDTraceType trace_types[] = {AGDMissing};
  struct BlendConfig blend_config = {};
  struct BlendParams blend_params = {};

  EXPECT_EQ(set_blend_params(n_traces, n_times, shottimes, channels,
                             trace_types, &blend_config, &blend_params),
            0);
  EXPECT_EQ(blend_params.n_traces, 1);
  EXPECT_EQ(blend_params.trace_length, 50);
  EXPECT_EQ(blend_params.blend_coords[0].channel_idx, -1);

  free_blend_params(&blend_params);
  free_blend_config(&blend_config);
}

TEST(Blend, ForwardAndAdjoint) {
  int n_patches = 2;
  int n_traces_per_patch[2] = {3, 4};
  int n_times_per_patch[2] = {3, 4};
  std::vector<std::vector<long int>> shottimes(2);
  std::vector<long int *> shottimes_ptr(n_patches, nullptr);
  std::vector<std::vector<int>> channels(2);
  std::vector<int *> channels_ptr(n_patches, nullptr);
  std::vector<std::vector<enum AGDTraceType>> trace_types(2);
  std::vector<enum AGDTraceType *> trace_types_ptr(n_patches, nullptr);
  std::vector<std::vector<AGD_TYPE>> data(2);
  std::vector<AGD_TYPE *> data_ptr(n_patches, nullptr);
  struct BlendConfig blend_config = {};
  std::vector<struct BlendParams> blend_params(n_patches);
  AGD_TYPE ***blended;

  shottimes[0].push_back(-5);
  shottimes[0].push_back(-3);
  shottimes[0].push_back(-1);
  shottimes[1].push_back(-10);
  shottimes[1].push_back(-4);
  shottimes[1].push_back(-2);
  shottimes[1].push_back(0);
  channels[0].resize(n_traces_per_patch[0]);
  std::fill(channels[0].begin(), channels[0].end(), 123);
  channels[1].resize(n_traces_per_patch[1]);
  std::fill(channels[1].begin(), channels[1].end(), 123);
  channels[1][0] = 456;
  trace_types[0].push_back(AGDMissing);
  trace_types[0].push_back(AGDLive);
  trace_types[0].push_back(AGDLive);
  trace_types[1].push_back(AGDLive);
  trace_types[1].push_back(AGDLive);
  trace_types[1].push_back(AGDBad);
  trace_types[1].push_back(AGDMissing);
  data[0].resize(n_traces_per_patch[0] * n_times_per_patch[0]);
  std::iota(data[0].begin(), data[0].end(), (AGD_TYPE)1.0);
  data[1].resize(n_traces_per_patch[1] * n_times_per_patch[1]);
  std::iota(data[1].begin(), data[1].end(), (AGD_TYPE)50.0);

  /* Patch 0: (3 traces, 3 times)
   * Channel -1:
   * 1X  2X  3X
   * Channel 0:
   *         4   5   6
   *                 7   8    9
   *
   * Patch 1: (4 traces, 4 times)
   * Channel 1:
   * 50  51  52  53
   * Channel 0:
   *     54  55  56  57
   *             58M 59M 60M 61M
   * Channel -1:
   * 62X 63X 64X 65X
   *
   * Blended Channel 0 before mute: (overwrite)
   *     54  55  58  59  60  61
   *
   * Blended Channel 0 before mute: (sum)
   *     54  59 119 129  68  70
   *
   * Blended Channel 0 after mute:
   *     54  59   0   0   0   0
   *
   * Channel 0 Adjoint:
   * Patch 0:
   *         59   0   0
   *                  0   0   0
   *
   * Patch 1:
   *     54  59   0   0   0   0
   *              0   0   0   0   0   0
   *
   */

  for (int i = 0; i < n_patches; i++) {
    shottimes_ptr[i] = &shottimes[i][0];
    channels_ptr[i] = &channels[i][0];
    trace_types_ptr[i] = &trace_types[i][0];
    data_ptr[i] = &data[i][0];
  }

  /* blend_config */
  EXPECT_EQ(set_blend_config(n_patches, n_traces_per_patch, n_times_per_patch,
                             shottimes_ptr.data(), channels_ptr.data(),
                             trace_types_ptr.data(), &blend_config),
            0);
  EXPECT_EQ(blend_config.n_channels, 2);
  EXPECT_EQ(blend_config.n_mute_coords, 1);
  EXPECT_EQ(blend_config.unique_channel_values[0], 123);
  EXPECT_EQ(blend_config.unique_channel_values[1], 456);
  EXPECT_EQ(blend_config.channels_intervals[0].n_intervals, 1);
  EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].start, -4);
  EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].stop, 2);
  EXPECT_EQ(blend_config.channels_intervals[1].n_intervals, 1);
  EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].start, -10);
  EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].stop, -6);
  EXPECT_EQ(blend_config.mute_coords[0].channel_idx, 0);
  EXPECT_EQ(blend_config.mute_coords[0].interval_idx, 0);
  EXPECT_EQ(blend_config.mute_coords[0].interval_start, 2);
  EXPECT_EQ(blend_config.mute_coords[0].n_times, 4);

  /* blend_params */
  for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    EXPECT_EQ(set_blend_params(
                  n_traces_per_patch[patch_idx], n_times_per_patch[patch_idx],
                  shottimes[patch_idx].data(), channels[patch_idx].data(),
                  trace_types[patch_idx].data(), &blend_config,
                  blend_params.data() + patch_idx),
              0);
    EXPECT_EQ(blend_params[patch_idx].n_traces, n_traces_per_patch[patch_idx]);
    EXPECT_EQ(blend_params[patch_idx].trace_length,
              n_times_per_patch[patch_idx]);
  }
  EXPECT_EQ(blend_params[0].blend_coords[0].channel_idx, -1);
  EXPECT_EQ(blend_params[0].blend_coords[1].channel_idx, 0);
  EXPECT_EQ(blend_params[0].blend_coords[1].interval_idx, 0);
  EXPECT_EQ(blend_params[0].blend_coords[1].interval_start, 1);
  EXPECT_EQ(blend_params[0].blend_coords[1].trace_start, 0);
  EXPECT_EQ(blend_params[0].blend_coords[1].n_times, 3);
  EXPECT_EQ(blend_params[0].blend_coords[2].channel_idx, 0);
  EXPECT_EQ(blend_params[0].blend_coords[2].interval_idx, 0);
  EXPECT_EQ(blend_params[0].blend_coords[2].interval_start, 3);
  EXPECT_EQ(blend_params[0].blend_coords[2].trace_start, 0);
  EXPECT_EQ(blend_params[0].blend_coords[2].n_times, 3);
  EXPECT_EQ(blend_params[1].blend_coords[0].channel_idx, 1);
  EXPECT_EQ(blend_params[1].blend_coords[0].interval_idx, 0);
  EXPECT_EQ(blend_params[1].blend_coords[0].interval_start, 0);
  EXPECT_EQ(blend_params[1].blend_coords[0].trace_start, 0);
  EXPECT_EQ(blend_params[1].blend_coords[0].n_times, 4);
  EXPECT_EQ(blend_params[1].blend_coords[1].channel_idx, 0);
  EXPECT_EQ(blend_params[1].blend_coords[1].interval_idx, 0);
  EXPECT_EQ(blend_params[1].blend_coords[1].interval_start, 0);
  EXPECT_EQ(blend_params[1].blend_coords[1].trace_start, 0);
  EXPECT_EQ(blend_params[1].blend_coords[1].n_times, 4);
  EXPECT_EQ(blend_params[1].blend_coords[2].channel_idx, 0);
  EXPECT_EQ(blend_params[1].blend_coords[2].interval_idx, 0);
  EXPECT_EQ(blend_params[1].blend_coords[2].interval_start, 2);
  EXPECT_EQ(blend_params[1].blend_coords[2].trace_start, 0);
  EXPECT_EQ(blend_params[1].blend_coords[2].n_times, 4);
  EXPECT_EQ(blend_params[1].blend_coords[3].channel_idx, -1);

  /* allocate blended */
  EXPECT_EQ(allocate_blended(&blend_config, &blended), 0);
  zero_blended(&blend_config, (AGD_TYPE *const *const *)blended);

  for (int channel_idx = 0; channel_idx < blend_config.n_channels;
       channel_idx++) {
    struct ChannelIntervals const *const channel_intervals =
        blend_config.channels_intervals + channel_idx;
    for (int interval_idx = 0; interval_idx < channel_intervals->n_intervals;
         interval_idx++) {
      struct Interval const *const interval =
          channel_intervals->intervals + interval_idx;
      for (int time_idx = 0; time_idx < interval->stop - interval->start;
           time_idx++) {
        EXPECT_FLOAT_EQ(blended[channel_idx][interval_idx][time_idx], AGD_ZERO);
      }
    }
  }

  /* blend overwrite */
  for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    blend_overwrite_forward(data[patch_idx].data(),
                            blend_params.data() + patch_idx,
                            (AGD_TYPE *const *const *)blended);
  }

  EXPECT_FLOAT_EQ(blended[0][0][0], (AGD_TYPE)54.0);
  EXPECT_FLOAT_EQ(blended[0][0][1], (AGD_TYPE)55.0);
  EXPECT_FLOAT_EQ(blended[0][0][2], (AGD_TYPE)58.0);
  EXPECT_FLOAT_EQ(blended[0][0][3], (AGD_TYPE)59.0);
  EXPECT_FLOAT_EQ(blended[0][0][4], (AGD_TYPE)60.0);
  EXPECT_FLOAT_EQ(blended[0][0][5], (AGD_TYPE)61.0);
  EXPECT_FLOAT_EQ(blended[1][0][0], (AGD_TYPE)50.0);
  EXPECT_FLOAT_EQ(blended[1][0][1], (AGD_TYPE)51.0);
  EXPECT_FLOAT_EQ(blended[1][0][2], (AGD_TYPE)52.0);
  EXPECT_FLOAT_EQ(blended[1][0][3], (AGD_TYPE)53.0);

  zero_blended(&blend_config, (AGD_TYPE *const *const *)blended);

  /* blend sum */
  for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    blend_sum_forward(data[patch_idx].data(), blend_params.data() + patch_idx,
                      (AGD_TYPE *const *const *)blended);
  }

  EXPECT_FLOAT_EQ(blended[0][0][0], (AGD_TYPE)54.0);
  EXPECT_FLOAT_EQ(blended[0][0][1], (AGD_TYPE)59.0);
  EXPECT_FLOAT_EQ(blended[0][0][2], (AGD_TYPE)119.0);
  EXPECT_FLOAT_EQ(blended[0][0][3], (AGD_TYPE)129.0);
  EXPECT_FLOAT_EQ(blended[0][0][4], (AGD_TYPE)68.0);
  EXPECT_FLOAT_EQ(blended[0][0][5], (AGD_TYPE)70.0);
  EXPECT_FLOAT_EQ(blended[1][0][0], (AGD_TYPE)50.0);
  EXPECT_FLOAT_EQ(blended[1][0][1], (AGD_TYPE)51.0);
  EXPECT_FLOAT_EQ(blended[1][0][2], (AGD_TYPE)52.0);
  EXPECT_FLOAT_EQ(blended[1][0][3], (AGD_TYPE)53.0);

  apply_mute(&blend_config, (AGD_TYPE *const *const *)blended);

  EXPECT_FLOAT_EQ(blended[0][0][0], (AGD_TYPE)54.0);
  EXPECT_FLOAT_EQ(blended[0][0][1], (AGD_TYPE)59.0);
  EXPECT_FLOAT_EQ(blended[0][0][2], (AGD_TYPE)0.0);
  EXPECT_FLOAT_EQ(blended[0][0][3], (AGD_TYPE)0.0);
  EXPECT_FLOAT_EQ(blended[0][0][4], (AGD_TYPE)0.0);
  EXPECT_FLOAT_EQ(blended[0][0][5], (AGD_TYPE)0.0);
  EXPECT_FLOAT_EQ(blended[1][0][0], (AGD_TYPE)50.0);
  EXPECT_FLOAT_EQ(blended[1][0][1], (AGD_TYPE)51.0);
  EXPECT_FLOAT_EQ(blended[1][0][2], (AGD_TYPE)52.0);
  EXPECT_FLOAT_EQ(blended[1][0][3], (AGD_TYPE)53.0);

  for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    blend_adjoint((AGD_TYPE const *const *const *)blended,
                  blend_params.data() + patch_idx, data[patch_idx].data());
  }

  for (int time_idx = 0; time_idx < 3; time_idx++) {
    EXPECT_FLOAT_EQ(data[0][0 * 3 + time_idx], AGD_ZERO);
  }
  EXPECT_FLOAT_EQ(data[0][1 * 3 + 0], (AGD_TYPE)59.0);
  for (int time_idx = 1; time_idx < 3; time_idx++) {
    EXPECT_FLOAT_EQ(data[0][1 * 3 + time_idx], AGD_ZERO);
  }
  for (int time_idx = 0; time_idx < 3; time_idx++) {
    EXPECT_FLOAT_EQ(data[0][2 * 3 + time_idx], AGD_ZERO);
  }
  EXPECT_FLOAT_EQ(data[1][0 * 4 + 0], (AGD_TYPE)50.0);
  EXPECT_FLOAT_EQ(data[1][0 * 4 + 1], (AGD_TYPE)51.0);
  EXPECT_FLOAT_EQ(data[1][0 * 4 + 2], (AGD_TYPE)52.0);
  EXPECT_FLOAT_EQ(data[1][0 * 4 + 3], (AGD_TYPE)53.0);
  EXPECT_FLOAT_EQ(data[1][1 * 4 + 0], (AGD_TYPE)54.0);
  EXPECT_FLOAT_EQ(data[1][1 * 4 + 1], (AGD_TYPE)59.0);
  for (int time_idx = 2; time_idx < 4; time_idx++) {
    EXPECT_FLOAT_EQ(data[1][1 * 4 + time_idx], AGD_ZERO);
  }
  for (int trace_idx = 2; trace_idx < 4; trace_idx++) {
    for (int time_idx = 0; time_idx < 4; time_idx++) {
      EXPECT_FLOAT_EQ(data[1][trace_idx * 4 + time_idx], AGD_ZERO);
    }
  }

  free_blended(&blend_config, &blended);
  for (auto &blend_param : blend_params) {
    free_blend_params(&blend_param);
  }
  free_blend_config(&blend_config);
}

TEST(Blend, MeanForward) {
  int n_patches = 2;
  int n_traces_per_patch[2] = {2, 2};
  int n_times_per_patch[2] = {10, 11};
  int taper_length = 3;
  AGD_TYPE weight[3] = {AGD_HALF * AGD_SQRT(AGD_TWO - AGD_SQRT(AGD_TWO)),
                        AGD_HALF * AGD_SQRT(AGD_TWO),
                        AGD_HALF * AGD_SQRT(AGD_TWO + AGD_SQRT(AGD_TWO))};
  std::vector<std::vector<long int>> shottimes(2);
  std::vector<long int *> shottimes_ptr(n_patches, nullptr);
  std::vector<std::vector<int>> channels(2);
  std::vector<int *> channels_ptr(n_patches, nullptr);
  std::vector<std::vector<enum AGDTraceType>> trace_types(2);
  std::vector<enum AGDTraceType *> trace_types_ptr(n_patches, nullptr);
  std::vector<std::vector<AGD_TYPE>> data(2);
  std::vector<AGD_TYPE *> data_ptr(n_patches, nullptr);
  struct BlendConfig blend_config = {};
  std::vector<struct BlendParams> blend_params(n_patches);
  AGD_TYPE ***blended;

  shottimes[0].push_back(-5);
  shottimes[0].push_back(-3);
  shottimes[1].push_back(-10);
  shottimes[1].push_back(-4);
  channels[0].resize(n_traces_per_patch[0]);
  std::fill(channels[0].begin(), channels[0].end(), 123);
  channels[1].resize(n_traces_per_patch[1]);
  std::fill(channels[1].begin(), channels[1].end(), 123);
  channels[1][0] = 456;
  trace_types[0].push_back(AGDLive);
  trace_types[0].push_back(AGDLive);
  trace_types[1].push_back(AGDLive);
  trace_types[1].push_back(AGDLive);
  data[0].resize(n_traces_per_patch[0] * n_times_per_patch[0]);
  std::iota(data[0].begin(), data[0].end(), (AGD_TYPE)1.0);
  data[1].resize(n_traces_per_patch[1] * n_times_per_patch[1]);
  std::iota(data[1].begin(), data[1].end(), (AGD_TYPE)50.0);

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

  std::vector<AGD_TYPE> expected_chan0_count{weight[0],
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
                             trace_types_ptr.data(), &blend_config),
            0);
  EXPECT_EQ(blend_config.n_channels, 2);
  EXPECT_EQ(blend_config.n_mute_coords, 0);
  EXPECT_EQ(blend_config.unique_channel_values[0], 123);
  EXPECT_EQ(blend_config.unique_channel_values[1], 456);
  EXPECT_EQ(blend_config.channels_intervals[0].n_intervals, 1);
  EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].start, -5);
  EXPECT_EQ(blend_config.channels_intervals[0].intervals[0].stop, 7);
  EXPECT_EQ(blend_config.channels_intervals[1].n_intervals, 1);
  EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].start, -10);
  EXPECT_EQ(blend_config.channels_intervals[1].intervals[0].stop, 1);

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
    EXPECT_FLOAT_EQ(blended[0][0][time_idx],
                    expected_chan0[time_idx] / expected_chan0_count[time_idx]);
  }
  for (int time_idx = 0; time_idx < 11; time_idx++) {
    EXPECT_FLOAT_EQ(blended[1][0][time_idx], data[1][time_idx]);
  }

  free_blended(&blend_config, &blended);
  for (auto &blend_param : blend_params) {
    free_blend_params(&blend_param);
  }
  free_blend_config(&blend_config);
}

TEST(Blend, DotTest) {
  srand(1);
  for (int rand_repeat = 0; rand_repeat < 5; rand_repeat++) {
    int n_patches = rand() % 5;
    std::vector<int> n_traces_per_patch(n_patches);
    std::vector<int> n_times_per_patch(n_patches);
    std::vector<std::vector<long int>> shottimes;
    std::vector<long int *> shottimes_ptr(n_patches, nullptr);
    std::vector<std::vector<int>> channels;
    std::vector<int *> channels_ptr(n_patches, nullptr);
    std::vector<std::vector<enum AGDTraceType>> trace_types;
    std::vector<enum AGDTraceType *> trace_types_ptr(n_patches, nullptr);
    std::vector<std::vector<AGD_TYPE>> x;
    std::vector<std::vector<AGD_TYPE>> xp;
    struct BlendConfig blend_config = {};
    std::vector<struct BlendParams> blend_params(n_patches);
    AGD_TYPE ***y;
    AGD_TYPE ***yp;

    int max_n_traces_per_patch = 64;
    std::generate(
        n_traces_per_patch.begin(), n_traces_per_patch.end(),
        [max_n_traces_per_patch]() { return rand() % max_n_traces_per_patch; });

    int max_n_times_per_patch = 64;
    std::generate(
        n_times_per_patch.begin(), n_times_per_patch.end(),
        [max_n_times_per_patch]() { return rand() % max_n_times_per_patch; });

    int shottime_max = std::inner_product(n_traces_per_patch.cbegin(),
                                          n_traces_per_patch.cend(),
                                          n_times_per_patch.cbegin(), 0) /
                       3;
    int channel_max = std::accumulate(n_traces_per_patch.cbegin(),
                                      n_traces_per_patch.cend(), 0) /
                      10;

    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      shottimes.push_back(std::vector<long int>(n_traces_per_patch[patch_idx]));
      std::generate(shottimes[patch_idx].begin(), shottimes[patch_idx].end(),
                    [shottime_max]() { return rand() % shottime_max; });

      channels.push_back(std::vector<int>(n_traces_per_patch[patch_idx]));
      std::generate(channels[patch_idx].begin(), channels[patch_idx].end(),
                    [channel_max]() { return rand() % channel_max; });

      trace_types.push_back(
          std::vector<enum AGDTraceType>(n_traces_per_patch[patch_idx]));
      std::generate(trace_types[patch_idx].begin(),
                    trace_types[patch_idx].end(),
                    []() { return (enum AGDTraceType)(rand() % 3); });

      x.push_back(std::vector<AGD_TYPE>(n_traces_per_patch[patch_idx] *
                                        n_times_per_patch[patch_idx]));
      std::generate(x[patch_idx].begin(), x[patch_idx].end(),
                    []() { return (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX; });

      xp.push_back(std::vector<AGD_TYPE>(n_traces_per_patch[patch_idx] *
                                         n_times_per_patch[patch_idx]));
      std::fill(xp[patch_idx].begin(), xp[patch_idx].end(), AGD_ZERO);
    }

    for (int i = 0; i < n_patches; i++) {
      shottimes_ptr[i] = &shottimes[i][0];
      channels_ptr[i] = &channels[i][0];
      trace_types_ptr[i] = &trace_types[i][0];
    }

    /* blend_config */
    EXPECT_EQ(set_blend_config(n_patches, n_traces_per_patch.data(),
                               n_times_per_patch.data(), shottimes_ptr.data(),
                               channels_ptr.data(), trace_types_ptr.data(),
                               &blend_config),
              0);

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
    EXPECT_EQ(allocate_blended(&blend_config, &y), 0);
    EXPECT_EQ(allocate_blended(&blend_config, &yp), 0);
    zero_blended(&blend_config, (AGD_TYPE *const *const *)y);

    for (int channel_idx = 0; channel_idx < blend_config.n_channels;
         channel_idx++) {
      struct ChannelIntervals const *const channel_intervals =
          blend_config.channels_intervals + channel_idx;
      for (int interval_idx = 0; interval_idx < channel_intervals->n_intervals;
           interval_idx++) {
        struct Interval const *const interval =
            channel_intervals->intervals + interval_idx;
        for (int time_idx = 0; time_idx < interval->stop - interval->start;
             time_idx++) {
          yp[channel_idx][interval_idx][time_idx] =
              (AGD_TYPE)rand() / (AGD_TYPE)RAND_MAX;
        }
      }
    }

    /* blend */
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      blend_sum_forward(x[patch_idx].data(), blend_params.data() + patch_idx,
                        (AGD_TYPE *const *const *)y);
    }

    apply_mute(&blend_config, (AGD_TYPE *const *const *)y);
    apply_mute(&blend_config, (AGD_TYPE *const *const *)yp);

    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      blend_adjoint((AGD_TYPE const *const *const *)yp,
                    blend_params.data() + patch_idx, xp[patch_idx].data());
    }

    AGD_TYPE sum_x = AGD_ZERO;
    for (int patch_idx = 0; patch_idx < n_patches; patch_idx++) {
      sum_x += std::inner_product(x[patch_idx].cbegin(), x[patch_idx].cend(),
                                  xp[patch_idx].begin(), AGD_ZERO);
    }

    AGD_TYPE sum_y = AGD_ZERO;
    for (int channel_idx = 0; channel_idx < blend_config.n_channels;
         channel_idx++) {
      struct ChannelIntervals const *const channel_intervals =
          blend_config.channels_intervals + channel_idx;
      for (int interval_idx = 0; interval_idx < channel_intervals->n_intervals;
           interval_idx++) {
        struct Interval const *const interval =
            channel_intervals->intervals + interval_idx;
        int n_times = interval->stop - interval->start;
        sum_y +=
            std::inner_product(y[channel_idx][interval_idx] + 0,
                               y[channel_idx][interval_idx] + n_times,
                               yp[channel_idx][interval_idx] + 0, AGD_ZERO);
      }
    }

#ifdef AGD_DOUBLE
    EXPECT_NEAR(sum_x, sum_y, 1e-12);
#else
    EXPECT_NEAR(sum_x, sum_y, 1e-4);
#endif

    free_blended(&blend_config, &y);
    free_blended(&blend_config, &yp);
    for (auto &blend_param : blend_params) {
      free_blend_params(&blend_param);
    }
    free_blend_config(&blend_config);
  }
}
