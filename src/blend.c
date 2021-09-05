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

static int interval_overlaps(struct Interval const *const interval,
                             struct Interval const *const other) {
  return (interval->start < other->stop) && (interval->stop > other->start);
}

static int interval_adjacent(struct Interval const *const interval,
                             struct Interval const *const other) {
  return (interval->start == other->stop) || (interval->stop == other->start);
}

#ifdef AGD_MPI
static int get_overlapping_interval_idx(
    struct Interval const *const interval,
    struct ChannelIntervals const *const channel_intervals,
    int const start_idx) {
  int interval_idx;

  for (interval_idx = start_idx; interval_idx < channel_intervals->n_intervals;
       interval_idx++) {
    struct Interval const *const existing_interval =
        channel_intervals->intervals + interval_idx;
    if (interval_overlaps(interval, existing_interval)) {
      return interval_idx;
    }
  }

  return -1;
}
#endif /* AGD_MPI */

static int get_overlapping_or_adjacent_interval_idx(
    struct Interval const *const interval,
    struct ChannelIntervals const *const channel_intervals,
    int const start_idx) {
  int interval_idx;

  for (interval_idx = start_idx; interval_idx < channel_intervals->n_intervals;
       interval_idx++) {
    struct Interval const *const existing_interval =
        channel_intervals->intervals + interval_idx;
    if (interval_overlaps(interval, existing_interval) ||
        interval_adjacent(interval, existing_interval)) {
      return interval_idx;
    }
  }

  return -1;
}

static void extend_interval(struct Interval *const interval,
                            struct Interval const *const other) {
  if (other->start < interval->start) {
    interval->start = other->start;
  }
  if (other->stop > interval->stop) {
    interval->stop = other->stop;
  }
}

static void delete_interval_idx(
    int const overlapping_idx,
    struct ChannelIntervals *const channel_intervals) {
  int interval_idx;

  /* Move up later intervals and then decrease the number of intervals.
   * This does not change the amount of allocated memory, but this
   * should not be a problem. */
  for (interval_idx = overlapping_idx;
       interval_idx < channel_intervals->n_intervals - 1; interval_idx++) {
    channel_intervals->intervals[interval_idx] =
        channel_intervals->intervals[interval_idx + 1];
  }

  channel_intervals->n_intervals--;
}

static int append_to_channel_intervals(
    struct Interval const *const interval,
    struct ChannelIntervals *const channel_intervals) {
  struct Interval *ptr;
  channel_intervals->n_intervals++;
  ptr = (struct Interval *)realloc(
      channel_intervals->intervals,
      (size_t)channel_intervals->n_intervals * sizeof(struct Interval));
  if (ptr == NULL) goto err;
  channel_intervals->intervals = ptr;
  channel_intervals->intervals[channel_intervals->n_intervals - 1] = *interval;

  return 0;
err:
  fprintf(stderr, "ERROR in append_to_channel_intervals\n");
  return 1;
}

static int add_interval(struct Interval const *const interval,
                        struct ChannelIntervals *const channel_intervals) {
  /* Check if overlaps with existing intervals
   * Begin by finding the first overlapping interval */
  int start_idx = 0;
  int const first_overlapping_idx = get_overlapping_or_adjacent_interval_idx(
      interval, channel_intervals, start_idx);

  if (first_overlapping_idx >= 0) {
    /* The new interval overlaps with at least one existing interval */
    struct Interval *const first_overlapping =
        channel_intervals->intervals + first_overlapping_idx;

    /* Extend (if needed) to include new interval */
    extend_interval(first_overlapping, interval);

    /* Now find any other intervals that overlap,
     * extend the first interval to include them,
     * and then delete them */
    start_idx = first_overlapping_idx + 1;
    while (1) {
      int overlapping_idx = get_overlapping_or_adjacent_interval_idx(
          interval, channel_intervals, start_idx);
      if (overlapping_idx < 0) break; /* No more overlapping found */
      extend_interval(first_overlapping,
                      channel_intervals->intervals + overlapping_idx);
      delete_interval_idx(overlapping_idx, channel_intervals);
      start_idx = overlapping_idx + 1;
    }
  } else {
    /* No overlapping interval found, so append as new interval */
    if (append_to_channel_intervals(interval, channel_intervals)) goto err;
  }

  return 0;
err:
  fprintf(stderr, "ERROR in add_interval\n");
  return 1;
}

static void free_channels_intervals(
    int const n_channels, struct ChannelIntervals **const channels_intervals) {
  if (*channels_intervals != NULL) {
    int channel_idx;
    for (channel_idx = 0; channel_idx < n_channels; channel_idx++) {
      free((*channels_intervals)[channel_idx].intervals);
      (*channels_intervals)[channel_idx].intervals = NULL;
    }
    free(*channels_intervals);
    *channels_intervals = NULL;
  }
}

/* Make a list of the unique channel values
 *
 * e.g., [[123, 456], [-789, 123]] -> [123, 456, -789]
 *
 * Channels with trace type AGDMissing will be ignored
 *
 * `unique_channel_values` should be a pointer to NULL.
 * It will be allocated during the function and so
 * will need to be freed by the caller.
 */
static int set_unique_channel_values(
    int const n_patches, int const *const n_traces,
    int const *const *const channels,
    enum AGDTraceType const *const *const trace_types,
    int *const n_unique_channels, int **const unique_channel_values) {
  int patch_idx;
  *n_unique_channels = 0;

  /* Find the unique channels in each patch */
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const patch_n_channels = n_traces[patch_idx];
    int const *const patch_channels = channels[patch_idx];
    enum AGDTraceType const *const patch_trace_types = trace_types[patch_idx];
    int channel_idx;

    for (channel_idx = 0; channel_idx < patch_n_channels; channel_idx++) {
      int const channel = patch_channels[channel_idx];
      int unique_idx;

      if (patch_trace_types[channel_idx] == AGDMissing) {
        continue;
      }

      for (unique_idx = 0; unique_idx < *n_unique_channels; unique_idx++) {
        if ((*unique_channel_values)[unique_idx] == channel) {
          /* Already in unique_channel_values list */
          break;
        }
      }

      if (unique_idx == *n_unique_channels) {
        /* A new unique channel value */
        int *ptr;
        (*n_unique_channels)++;
        ptr = (int *)realloc(*unique_channel_values,
                             (size_t)*n_unique_channels * sizeof(int));
        if (ptr == NULL) goto err;
        *unique_channel_values = ptr;
        (*unique_channel_values)[*n_unique_channels - 1] = channel;
      }
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_unique_channel_values\n");
  return 1;
}

/* Convert channel numbers to their index in unique_channel_values
 *
 * The index of channels with trace type AGDMissing, and channels
 * that are not in `unique_channel_values`, will be -1.
 *
 */
static int get_channel_idx(int const channel,
                           enum AGDTraceType const trace_type,
                           int const n_unique_channels,
                           int const *const unique_channel_values) {
  int unique_idx;

  if (trace_type == AGDMissing) return -1;

  for (unique_idx = 0; unique_idx < n_unique_channels; unique_idx++) {
    if (unique_channel_values[unique_idx] == channel) {
      return unique_idx;
    }
  }

  return -1;
}

static int get_n_mute_coords(
    int const n_channels,
    struct ChannelIntervals const *channels_mute_intervals) {
  int n_mute_coords = 0;
  int channel_idx;
  for (channel_idx = 0; channel_idx < n_channels; channel_idx++) {
    struct ChannelIntervals const *const channel_mute_intervals =
        channels_mute_intervals + channel_idx;
    n_mute_coords += channel_mute_intervals->n_intervals;
  }
  return n_mute_coords;
}

static void free_blend_config(struct BlendConfig *const blend_config) {
  free_channels_intervals(blend_config->n_channels,
                          &(blend_config->channels_intervals));
  free(blend_config->unique_channel_values);
  blend_config->unique_channel_values = NULL;
  free(blend_config->mute_coords);
  blend_config->mute_coords = NULL;
#ifdef AGD_MPI
  if (blend_config->ranks_coords != NULL) {
    int rank_idx;
    for (rank_idx = 0; rank_idx < blend_config->comm_size; rank_idx++) {
      free(blend_config->ranks_coords[rank_idx]);
      blend_config->ranks_coords[rank_idx] = NULL;
    }
    free(blend_config->ranks_coords);
    blend_config->ranks_coords = NULL;
  }
  if (blend_config->send_buffers != NULL) {
    int rank_idx;
    for (rank_idx = 0; rank_idx < blend_config->comm_size; rank_idx++) {
      free(blend_config->send_buffers[rank_idx]);
      blend_config->send_buffers[rank_idx] = NULL;
    }
    free(blend_config->send_buffers);
    blend_config->send_buffers = NULL;
  }
  if (blend_config->receive_buffers != NULL) {
    int rank_idx;
    for (rank_idx = 0; rank_idx < blend_config->comm_size; rank_idx++) {
      free(blend_config->receive_buffers[rank_idx]);
      blend_config->receive_buffers[rank_idx] = NULL;
    }
    free(blend_config->receive_buffers);
    blend_config->receive_buffers = NULL;
  }
  free(blend_config->requests);
  blend_config->requests = NULL;
  free(blend_config->n_overlaps);
  blend_config->n_overlaps = NULL;
#endif /* AGD_MPI */
  blend_config->n_channels = 0;
  blend_config->n_mute_coords = 0;
}

static int set_channels_intervals(
    int const n_patches, int const n_channels, int const *const n_traces,
    int const *const n_times, long int const *const *const shottimes,
    int const *const *const channels,
    enum AGDTraceType const *const *const trace_types,
    int const *const unique_channel_values,
    struct ChannelIntervals *const channels_intervals,
    struct ChannelIntervals *const channels_mute_intervals) {
  int patch_idx;
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int const patch_n_traces = n_traces[patch_idx];
    int const patch_n_times = n_times[patch_idx];
    long int const *const patch_shottimes = shottimes[patch_idx];
    int const *const patch_channels = channels[patch_idx];
    enum AGDTraceType const *const patch_trace_types = trace_types[patch_idx];
    int trace_idx;
    for (trace_idx = 0; trace_idx < patch_n_traces; trace_idx++) {
      long int const shottime = patch_shottimes[trace_idx];
      enum AGDTraceType const trace_type = patch_trace_types[trace_idx];
      int const channel = get_channel_idx(patch_channels[trace_idx], trace_type,
                                          n_channels, unique_channel_values);
      struct Interval interval;
      interval.start = shottime;
      interval.stop = shottime + patch_n_times;
      if (channel < 0) continue; /* Missing */
      if (add_interval(&interval, channels_intervals + channel)) goto err;
      if (trace_type == AGDBad) {
        if (add_interval(&interval, channels_mute_intervals + channel))
          goto err;
      }
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_channels_intervals\n");
  return 1;
}

#ifdef AGD_MPI

/* Sort so that channel_idxs (slow) and interval_idxs (fast) are ascending */
static int rank_coords_sort(void const *const a, void const *const b) {
  struct BlendCoords const *const a2 = (struct BlendCoords const *)a;
  struct BlendCoords const *const b2 = (struct BlendCoords const *)b;
  if (a2->channel_idx < b2->channel_idx) {
    return -1;
  } else if (a2->channel_idx > b2->channel_idx) {
    return +1;
  }

  return a2->interval_idx - b2->interval_idx;
}

/* Setup so that can swap overlapping intervals of blended data.
 * This enables the blended data to be split over processes.
 *
 * Each rank broadcasts information about the intervals covered
 * by its data. Each rank receiving this information then determines
 * the overlaps with its own intervals. */
static int set_swap_intervals(
    int const n_channels,
    struct ChannelIntervals const *const channels_intervals,
    int const *const unique_channel_values, MPI_Comm comm,
    struct ChannelIntervals *const channels_mute_intervals,
    int **const n_overlaps, int *const comm_size, MPI_Request **const requests,
    AGD_TYPE ***const send_buffers, AGD_TYPE ***const receive_buffers,
    struct BlendCoords ***const ranks_coords) {
  int comm_rank;
  int rank_idx;
  if (MPI_Comm_rank(comm, &comm_rank)) goto err;
  if (MPI_Comm_size(comm, comm_size)) goto err;

  /* Allocate space to store information about overlapping intervals */
  *ranks_coords = (struct BlendCoords **)calloc((size_t)*comm_size,
                                                sizeof(struct BlendCoords *));
  if (*ranks_coords == NULL) goto err;
  *n_overlaps = (int *)calloc((size_t)*comm_size, sizeof(int));
  if (*n_overlaps == NULL) goto err;
  *send_buffers = (AGD_TYPE **)calloc((size_t)*comm_size, sizeof(AGD_TYPE *));
  if (*send_buffers == NULL) goto err;
  *receive_buffers =
      (AGD_TYPE **)calloc((size_t)*comm_size, sizeof(AGD_TYPE *));
  if (*receive_buffers == NULL) goto err;

  for (rank_idx = 0; rank_idx < *comm_size; rank_idx++) {
    int rank_n_channels;
    int *rank_unique_channel_values = NULL;
    size_t buffer_size = 0;
    int channel_idx;

    /* n_channels on rank */
    if (comm_rank == rank_idx) {
      rank_n_channels = n_channels;
    }
    if (MPI_Bcast(&rank_n_channels, 1, MPI_INT, rank_idx, comm)) goto err;

    /* unique_channel_values on rank */
    rank_unique_channel_values =
        (int *)malloc((size_t)rank_n_channels * sizeof(int));
    if (rank_unique_channel_values == NULL) goto err;
    if (comm_rank == rank_idx) {
      memcpy(rank_unique_channel_values, unique_channel_values,
             (size_t)n_channels * sizeof(int));
    }
    if (MPI_Bcast(rank_unique_channel_values, rank_n_channels, MPI_INT,
                  rank_idx, comm))
      goto err;

    for (channel_idx = 0; channel_idx < rank_n_channels; channel_idx++) {
      struct ChannelIntervals rank_intervals = ZERO_INIT;
      struct ChannelIntervals rank_mute_intervals = ZERO_INIT;

      /* n_intervals of rank's channel */
      if (comm_rank == rank_idx) {
        rank_intervals.n_intervals =
            channels_intervals[channel_idx].n_intervals;
        rank_mute_intervals.n_intervals =
            channels_mute_intervals[channel_idx].n_intervals;
      }
      if (MPI_Bcast(&(rank_intervals.n_intervals), 1, MPI_INT, rank_idx, comm))
        goto err;
      if (MPI_Bcast(&(rank_mute_intervals.n_intervals), 1, MPI_INT, rank_idx,
                    comm))
        goto err;

      /* intervals of rank's channel */
      rank_intervals.intervals = (struct Interval *)malloc(
          (size_t)rank_intervals.n_intervals * sizeof(struct Interval));
      if (rank_intervals.intervals == NULL) goto err;
      rank_mute_intervals.intervals = (struct Interval *)malloc(
          (size_t)rank_mute_intervals.n_intervals * sizeof(struct Interval));
      if (rank_mute_intervals.intervals == NULL) goto err;
      if (comm_rank == rank_idx) {
        memcpy(rank_intervals.intervals,
               channels_intervals[channel_idx].intervals,
               (size_t)rank_intervals.n_intervals * sizeof(struct Interval));
        memcpy(
            rank_mute_intervals.intervals,
            channels_mute_intervals[channel_idx].intervals,
            (size_t)rank_mute_intervals.n_intervals * sizeof(struct Interval));
      }
      if (MPI_Bcast(rank_intervals.intervals, 2 * rank_intervals.n_intervals,
                    MPI_LONG, rank_idx, comm))
        goto err;
      if (MPI_Bcast(rank_mute_intervals.intervals,
                    2 * rank_mute_intervals.n_intervals, MPI_LONG, rank_idx,
                    comm))
        goto err;

      if (comm_rank != rank_idx) {
        int local_channel_idx;

        /* Find local_channel_idx, if it exists */
        for (local_channel_idx = 0; local_channel_idx < n_channels;
             local_channel_idx++) {
          if (unique_channel_values[local_channel_idx] ==
              rank_unique_channel_values[channel_idx])
            break;
        }

        /* Loop over rank's channel intervals */
        if (local_channel_idx < n_channels) {
          struct ChannelIntervals const *const channel_intervals =
              channels_intervals + local_channel_idx;
          struct ChannelIntervals *const channel_mute_intervals =
              channels_mute_intervals + local_channel_idx;
          int interval_idx;
          for (interval_idx = 0; interval_idx < rank_intervals.n_intervals;
               interval_idx++) {
            /* If an interval overlaps with a local interval, add it to
             * ranks_coords */
            int local_interval_idx = get_overlapping_interval_idx(
                rank_intervals.intervals + interval_idx, channel_intervals, 0);
            while (local_interval_idx >= 0) {
              struct Interval const *const local_interval =
                  channel_intervals->intervals + local_interval_idx;
              struct Interval const *const rank_interval =
                  rank_intervals.intervals + interval_idx;
              struct BlendCoords *ptr = (struct BlendCoords *)realloc(
                  (*ranks_coords)[rank_idx],
                  (size_t)((*n_overlaps)[rank_idx] + 1) *
                      sizeof(struct BlendCoords));
              long int actual_start =
                  rank_interval->start < local_interval->start
                      ? local_interval->start
                      : rank_interval->start;
              long int actual_stop = rank_interval->stop > local_interval->stop
                                         ? local_interval->stop
                                         : rank_interval->stop;
              if (ptr == NULL) goto err;
              (*ranks_coords)[rank_idx] = ptr;
              (*n_overlaps)[rank_idx]++;
              ptr = (*ranks_coords)[rank_idx] + (*n_overlaps)[rank_idx] - 1;
              ptr->channel_idx = local_channel_idx;
              ptr->interval_idx = local_interval_idx;
              ptr->interval_start =
                  (size_t)(actual_start - local_interval->start);
              ptr->trace_start = (int)(actual_start - rank_interval->start);
              ptr->n_times = (int)(actual_stop - actual_start);
              buffer_size += (size_t)(ptr->n_times);
              local_interval_idx = get_overlapping_interval_idx(
                  rank_intervals.intervals + interval_idx, channel_intervals,
                  local_interval_idx + 1);
            }
          }

          /* Add rank's mutes to local mutes */
          for (interval_idx = 0; interval_idx < rank_mute_intervals.n_intervals;
               interval_idx++) {
            int local_interval_idx = get_overlapping_interval_idx(
                rank_mute_intervals.intervals + interval_idx, channel_intervals,
                0);
            /* Do not need to loop over get_overlapping...
               as we add the whole interval. Intersection with multiple
               intervals will be handled in set_mute_coords */
            if (local_interval_idx >= 0) {
              struct Interval const *const rank_interval =
                  rank_mute_intervals.intervals + interval_idx;
              if (add_interval(rank_interval, channel_mute_intervals)) goto err;
            }
          }
        }
      }

      /* Free temporary memory */
      free(rank_intervals.intervals);
      rank_intervals.intervals = NULL;
      free(rank_mute_intervals.intervals);
      rank_mute_intervals.intervals = NULL;
    }

    /* Swap the order of coords so can send/recv in the same order.
     * The order of channels and intervals on each rank may be different,
     * so if we did not do this then the order in which the samples are
     * received may be different from the order in which they are sent,
     * so we would need to store different orderings for packing and
     * unpacking.
     * We sort it so the order of the higher rank is always used. */
    if (comm_rank > rank_idx) {
      qsort((*ranks_coords)[rank_idx], (size_t)(*n_overlaps)[rank_idx],
            sizeof(struct BlendCoords), rank_coords_sort);
    }

    /* Allocate buffers */
    if (buffer_size > 0) {
      if (buffer_size >= 2147483647 / sizeof(AGD_TYPE)) {
        printf("WARNING from rank %d: %ld bytes will be swapped with rank %d\n",
               comm_rank, buffer_size * sizeof(AGD_TYPE), rank_idx);
      }
      if (buffer_size >= INT_MAX) {
        fprintf(stderr,
                "ERROR: There is too much blended overlap between ranks\n");
        goto err;
      }
      (*send_buffers)[rank_idx] =
          (AGD_TYPE *)malloc(buffer_size * sizeof(AGD_TYPE));
      if ((*send_buffers)[rank_idx] == NULL) goto err;
      (*receive_buffers)[rank_idx] =
          (AGD_TYPE *)malloc(buffer_size * sizeof(AGD_TYPE));
      if ((*receive_buffers)[rank_idx] == NULL) goto err;
    }

    /* Free temporary memory */
    free(rank_unique_channel_values);
    rank_unique_channel_values = NULL;
  }

  /* Allocate requests */
  *requests =
      (MPI_Request *)malloc((size_t)(2 * *comm_size) * sizeof(MPI_Request));
  if (*requests == NULL) goto err;

  return 0;
err:
  fprintf(stderr, "ERROR in set_swap_intervals\n");
  return 1;
}
#endif /* AGD_MPI */

static void set_mute_coords(
    int const n_channels,
    struct ChannelIntervals const *const channels_intervals,
    struct ChannelIntervals const *const channels_mute_intervals,
    struct BlendCoords *const mute_coords) {
  int channel_idx;
  int mute_coord_idx = 0;
  for (channel_idx = 0; channel_idx < n_channels; channel_idx++) {
    struct ChannelIntervals const *const channel_mute_intervals =
        channels_mute_intervals + channel_idx;
    int mute_interval_idx;
    for (mute_interval_idx = 0;
         mute_interval_idx < channel_mute_intervals->n_intervals;
         mute_interval_idx++) {
      struct Interval const *const mute_interval =
          channel_mute_intervals->intervals + mute_interval_idx;
      struct ChannelIntervals const *const channel_intervals =
          channels_intervals + channel_idx;
      int interval_idx;
      for (interval_idx = 0; interval_idx < channel_intervals->n_intervals;
           interval_idx++) {
        struct Interval const *const channel_interval =
            channel_intervals->intervals + interval_idx;
        if (interval_overlaps(mute_interval, channel_interval)) {
          struct BlendCoords *const mute_coord = mute_coords + mute_coord_idx;
          long int actual_start = mute_interval->start < channel_interval->start
                                      ? channel_interval->start
                                      : mute_interval->start;
          long int actual_stop = mute_interval->stop > channel_interval->stop
                                     ? channel_interval->stop
                                     : mute_interval->stop;
          mute_coord->channel_idx = channel_idx;
          mute_coord->interval_idx = interval_idx;
          mute_coord->interval_start =
              (size_t)(actual_start - channel_interval->start);
          mute_coord->trace_start = 0;
          mute_coord->n_times = (int)(actual_stop - actual_start);
          mute_coord_idx++;
        }
      }
    }
  }
}

static int set_blend_config(int const n_patches, int const *const n_traces,
                            int const *const n_times,
                            long int const *const *const shottimes,
                            int const *const *const channels,
                            enum AGDTraceType const *const *const trace_types,
#ifdef AGD_MPI
                            MPI_Comm comm,
#endif /* AGD_MPI */
                            struct BlendConfig *blend_config) {
  struct ChannelIntervals *channels_mute_intervals = NULL;

  /* Get the unique channel values */
  if (set_unique_channel_values(n_patches, n_traces, channels, trace_types,
                                &(blend_config->n_channels),
                                &(blend_config->unique_channel_values)))
    goto err;

  /* Allocate and initialize channel intervals */
  blend_config->channels_intervals = (struct ChannelIntervals *)calloc(
      (size_t)blend_config->n_channels, sizeof(struct ChannelIntervals));
  if (blend_config->channels_intervals == NULL) goto err;
  channels_mute_intervals = (struct ChannelIntervals *)calloc(
      (size_t)blend_config->n_channels, sizeof(struct ChannelIntervals));
  if (channels_mute_intervals == NULL) goto err;

  /* Set the intervals for each channel */
  if (set_channels_intervals(
          n_patches, blend_config->n_channels, n_traces, n_times, shottimes,
          channels, trace_types, blend_config->unique_channel_values,
          blend_config->channels_intervals, channels_mute_intervals))
    goto err;

#ifdef AGD_MPI
  /* Setup to swap overlapping blended intervals */
  if (set_swap_intervals(
          blend_config->n_channels, blend_config->channels_intervals,
          blend_config->unique_channel_values, comm, channels_mute_intervals,
          &(blend_config->n_overlaps), &(blend_config->comm_size),
          &(blend_config->requests), &(blend_config->send_buffers),
          &(blend_config->receive_buffers), &(blend_config->ranks_coords)))
    goto err;
  blend_config->comm = comm;
#endif /* AGD_MPI */

  /* Get the number of mute coords */
  blend_config->n_mute_coords =
      get_n_mute_coords(blend_config->n_channels, channels_mute_intervals);

  /* Allocate mute_coords */
  blend_config->mute_coords = (struct BlendCoords *)malloc(
      (size_t)blend_config->n_mute_coords * sizeof(struct BlendCoords));

  /* Convert the channel_mute_intervals into BlendCoords */
  set_mute_coords(blend_config->n_channels, blend_config->channels_intervals,
                  channels_mute_intervals, blend_config->mute_coords);

  free_channels_intervals(blend_config->n_channels, &channels_mute_intervals);
  return 0;
err:
  fprintf(stderr, "ERROR in set_blend_config\n");
  free_channels_intervals(blend_config->n_channels, &channels_mute_intervals);
  return 1;
}

static void set_blend_coords(
    int const n_times, long int const shottime, int const channel,
    enum AGDTraceType const trace_type,
    struct ChannelIntervals const *const channel_intervals,
    struct BlendCoords *const blend_coords) {
  struct Interval interval;
  interval.start = shottime;
  interval.stop = shottime + n_times;

  if ((trace_type != AGDMissing) && (channel >= 0)) {
    int interval_idx;
    blend_coords->channel_idx = channel;
    for (interval_idx = 0; interval_idx < channel_intervals->n_intervals;
         interval_idx++) {
      struct Interval const *const channel_interval =
          channel_intervals->intervals + interval_idx;
      if (interval_overlaps(&interval, channel_interval)) {
        long int actual_start = interval.start < channel_interval->start
                                    ? channel_interval->start
                                    : interval.start;
        long int actual_stop = interval.stop > channel_interval->stop
                                   ? channel_interval->stop
                                   : interval.stop;
        blend_coords->interval_idx = interval_idx;
        blend_coords->interval_start =
            (size_t)(actual_start - channel_interval->start);
        blend_coords->trace_start = (int)(actual_start - interval.start);
        blend_coords->n_times = (int)(actual_stop - actual_start);
        return;
      }
    }
  }
  blend_coords->channel_idx = -1;
}

static int set_blend_params(int const n_traces, int const n_times,
                            long int const *const shottimes,
                            int const *const channels,
                            enum AGDTraceType const *const trace_types,
                            struct BlendConfig const *const blend_config,
                            struct BlendParams *const blend_params) {
  int trace_idx;

  blend_params->n_traces = n_traces;
  blend_params->trace_length = n_times;

  /* Allocate blend_coords */
  blend_params->blend_coords = (struct BlendCoords *)malloc(
      (size_t)n_traces * sizeof(struct BlendCoords));
  if (blend_params->blend_coords == NULL) goto err;

  /* Set the coordinates of each trace in the blended data */
  for (trace_idx = 0; trace_idx < n_traces; trace_idx++) {
    long int const shottime = shottimes[trace_idx];
    enum AGDTraceType const trace_type = trace_types[trace_idx];
    int const channel = get_channel_idx(channels[trace_idx], trace_type,
                                        blend_config->n_channels,
                                        blend_config->unique_channel_values);
    set_blend_coords(n_times, shottime, channel, trace_type,
                     blend_config->channels_intervals + channel,
                     blend_params->blend_coords + trace_idx);
  }

  return 0;
err:
  fprintf(stderr, "ERROR in set_blend_params\n");
  return 1;
}

static void free_blend_params(struct BlendParams *const blend_params) {
  free(blend_params->blend_coords);
  blend_params->blend_coords = NULL;
}

static void free_blend_params_array(int const n_elements,
                                    struct BlendParams **const blend_params) {
  if (*blend_params != NULL) {
    int element_idx;
    for (element_idx = 0; element_idx < n_elements; element_idx++) {
      free_blend_params((*blend_params) + element_idx);
    }
    free(*blend_params);
    *blend_params = NULL;
  }
}

static void zero_blended(struct BlendConfig const *const blend_config,
                         AGD_TYPE *const *const *const blended) {
  int channel_idx;
  for (channel_idx = 0; channel_idx < blend_config->n_channels; channel_idx++) {
    struct ChannelIntervals const *const channel_intervals =
        blend_config->channels_intervals + channel_idx;
    int interval_idx;
    for (interval_idx = 0; interval_idx < channel_intervals->n_intervals;
         interval_idx++) {
      struct Interval const *const interval =
          channel_intervals->intervals + interval_idx;
      size_t n_times = (size_t)(interval->stop - interval->start);
      memset(blended[channel_idx][interval_idx], 0, n_times * sizeof(AGD_TYPE));
    }
  }
}

static int allocate_blended(struct BlendConfig const *const blend_config,
                            AGD_TYPE ****const blended) {
  int const n_channels = blend_config->n_channels;
  int channel_idx;

  /* Allocate dataset */
  *blended = (AGD_TYPE ***)calloc((size_t)n_channels, sizeof(AGD_TYPE **));
  if (*blended == NULL) goto err;

  /* Allocate channels */
  for (channel_idx = 0; channel_idx < n_channels; channel_idx++) {
    struct ChannelIntervals const *const channel_intervals =
        blend_config->channels_intervals + channel_idx;
    int const n_intervals = channel_intervals->n_intervals;
    int interval_idx;
    (*blended)[channel_idx] =
        (AGD_TYPE **)calloc((size_t)n_intervals, sizeof(AGD_TYPE *));
    if ((*blended)[channel_idx] == NULL) goto err;

    /* Allocate intervals */
    for (interval_idx = 0; interval_idx < n_intervals; interval_idx++) {
      struct Interval const *const interval =
          channel_intervals->intervals + interval_idx;
      (*blended)[channel_idx][interval_idx] = (AGD_TYPE *)malloc(
          (size_t)(interval->stop - interval->start) * sizeof(AGD_TYPE));
      if ((*blended)[channel_idx][interval_idx] == NULL) goto err;
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in allocate_blended\n");
  return 1;
}

static void free_blended(struct BlendConfig const *const blend_config,
                         AGD_TYPE ****const blended) {
  int n_channels = blend_config->n_channels;
  int channel_idx;

  if (*blended != NULL) {
    /* Free channels */
    for (channel_idx = 0; channel_idx < n_channels; channel_idx++) {
      struct ChannelIntervals const *const channel_intervals =
          blend_config->channels_intervals + channel_idx;
      int n_intervals = channel_intervals->n_intervals;
      int interval_idx;

      if ((*blended)[channel_idx] != NULL) {
        /* Free intervals */
        for (interval_idx = 0; interval_idx < n_intervals; interval_idx++) {
          free((*blended)[channel_idx][interval_idx]);
          (*blended)[channel_idx][interval_idx] = NULL;
        }

        free((*blended)[channel_idx]);
        (*blended)[channel_idx] = NULL;
      }
    }

    /* Free dataset */
    free(*blended);
    *blended = NULL;
  }
}

static void apply_mute(struct BlendConfig const *const blend_config,
                       AGD_TYPE *const *const *const blended) {
  int mute_idx;
  for (mute_idx = 0; mute_idx < blend_config->n_mute_coords; mute_idx++) {
    struct BlendCoords const *const mute_coords =
        blend_config->mute_coords + mute_idx;
    int const channel_idx = mute_coords->channel_idx;
    int const interval_idx = mute_coords->interval_idx;
    size_t const interval_start = mute_coords->interval_start;
    int const n_times = mute_coords->n_times;
    memset(blended[channel_idx][interval_idx] + interval_start, 0,
           (size_t)n_times * sizeof(AGD_TYPE));
  }
}

#ifdef AGD_MPI
static int blend_swap(struct BlendConfig const *const blend_config,
                      AGD_TYPE const *const *const *const blended) {
  MPI_Comm comm = blend_config->comm;
  int rank_idx;
  for (rank_idx = 0; rank_idx < blend_config->comm_size; rank_idx++) {
    if (blend_config->n_overlaps[rank_idx] > 0) {
      struct BlendCoords const *const rank_coords =
          blend_config->ranks_coords[rank_idx];
      AGD_TYPE *const send_buffer = blend_config->send_buffers[rank_idx];
      AGD_TYPE *const receive_buffer = blend_config->receive_buffers[rank_idx];
      int buffer_idx = 0;
      int overlap_idx;

      /* Copy to send buffer */
      for (overlap_idx = 0; overlap_idx < blend_config->n_overlaps[rank_idx];
           overlap_idx++) {
        struct BlendCoords const *const overlap_coords =
            rank_coords + overlap_idx;
        int const channel_idx = overlap_coords->channel_idx;
        int const interval_idx = overlap_coords->interval_idx;
        size_t const interval_start = overlap_coords->interval_start;
        int const n_times = overlap_coords->n_times;
        AGD_TYPE const *const blended_interval_start =
            blended[channel_idx][interval_idx] + interval_start;
        AGD_TYPE *const buffer_start = send_buffer + buffer_idx;
        memcpy(buffer_start, blended_interval_start,
               (size_t)n_times * sizeof(AGD_TYPE));
        buffer_idx += n_times;
      }

      /* Initiate send and receive */
      if (MPI_Isend(send_buffer, buffer_idx, AGD_MPI_TYPE, rank_idx, 0, comm,
                    blend_config->requests + rank_idx))
        goto err;
      if (MPI_Irecv(
              receive_buffer, buffer_idx, AGD_MPI_TYPE, rank_idx, MPI_ANY_TAG,
              comm,
              blend_config->requests + blend_config->comm_size + rank_idx))
        goto err;
    } else {
      blend_config->requests[rank_idx] = MPI_REQUEST_NULL;
      blend_config->requests[blend_config->comm_size + rank_idx] =
          MPI_REQUEST_NULL;
    }
  }

  /* Wait until sends and receives have finished */
  if (MPI_Waitall(2 * blend_config->comm_size, blend_config->requests,
                  MPI_STATUSES_IGNORE))
    goto err;

  return 0;
err:
  fprintf(stderr, "ERROR in blend_swap\n");
  return 1;
}
#endif /* AGD_MPI */

static void blend_sum_forward(AGD_TYPE const *const data,
                              struct BlendParams const *const blend_params,
                              AGD_TYPE *const *const *const blended) {
  int const trace_length = blend_params->trace_length;
  int trace_idx;
  for (trace_idx = 0; trace_idx < blend_params->n_traces; trace_idx++) {
    struct BlendCoords const *const trace_blend_coords =
        blend_params->blend_coords + trace_idx;
    int const channel_idx = trace_blend_coords->channel_idx;
    if (channel_idx >= 0) {
      int const interval_idx = trace_blend_coords->interval_idx;
      size_t const interval_start = trace_blend_coords->interval_start;
      int const trace_start = trace_blend_coords->trace_start;
      int const n_times = trace_blend_coords->n_times;
      int time_idx;
      AGD_TYPE *const blended_interval_start =
          blended[channel_idx][interval_idx] + interval_start;
      AGD_TYPE const *const data_start =
          data + (size_t)trace_idx * (size_t)trace_length + (size_t)trace_start;
      for (time_idx = 0; time_idx < n_times; time_idx++) {
        blended_interval_start[time_idx] += data_start[time_idx];
      }
    }
  }
}

#ifdef AGD_MPI
static int blend_sum_forward_mpi(struct BlendConfig const *const blend_config,
                                 AGD_TYPE *const *const *const blended) {
  int rank_idx;

  if (blend_swap(blend_config, (AGD_TYPE const *const *const *)blended))
    goto err;

  for (rank_idx = 0; rank_idx < blend_config->comm_size; rank_idx++) {
    if (blend_config->n_overlaps[rank_idx] > 0) {
      struct BlendCoords const *const rank_coords =
          blend_config->ranks_coords[rank_idx];
      AGD_TYPE const *const receive_buffer =
          blend_config->receive_buffers[rank_idx];
      int buffer_idx = 0;
      int overlap_idx;

      /* Sum receive buffer into blended */
      for (overlap_idx = 0; overlap_idx < blend_config->n_overlaps[rank_idx];
           overlap_idx++) {
        struct BlendCoords const *const overlap_coords =
            rank_coords + overlap_idx;
        int const channel_idx = overlap_coords->channel_idx;
        int const interval_idx = overlap_coords->interval_idx;
        size_t const interval_start = overlap_coords->interval_start;
        int const n_times = overlap_coords->n_times;
        int time_idx;
        AGD_TYPE *const blended_interval_start =
            blended[channel_idx][interval_idx] + interval_start;
        AGD_TYPE const *const buffer_start = receive_buffer + buffer_idx;
        for (time_idx = 0; time_idx < n_times; time_idx++) {
          blended_interval_start[time_idx] += buffer_start[time_idx];
        }
        buffer_idx += n_times;
      }
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in blend_sum_forward_mpi\n");
  return 1;
}
#endif /* AGD_MPI */

static void blend_overwrite_forward(
    AGD_TYPE const *const data, struct BlendParams const *const blend_params,
    AGD_TYPE *const *const *const blended) {
  int const trace_length = blend_params->trace_length;
  int trace_idx;
  for (trace_idx = 0; trace_idx < blend_params->n_traces; trace_idx++) {
    struct BlendCoords const *const trace_blend_coords =
        blend_params->blend_coords + trace_idx;
    int const channel_idx = trace_blend_coords->channel_idx;
    if (channel_idx >= 0) {
      int const interval_idx = trace_blend_coords->interval_idx;
      size_t const interval_start = trace_blend_coords->interval_start;
      int const trace_start = trace_blend_coords->trace_start;
      int const n_times = trace_blend_coords->n_times;
      memcpy(
          blended[channel_idx][interval_idx] + interval_start,
          data + (size_t)trace_idx * (size_t)trace_length + (size_t)trace_start,
          (size_t)n_times * sizeof(AGD_TYPE));
    }
  }
}

#ifdef AGD_MPI
static int blend_overwrite_forward_mpi(
    struct BlendConfig const *const blend_config,
    AGD_TYPE *const *const *const blended) {
  int comm_rank;
  int rank_idx;

  if (MPI_Comm_rank(blend_config->comm, &comm_rank)) goto err;

  if (blend_swap(blend_config, (AGD_TYPE const *const *const *)blended))
    goto err;

  /* To ensure that the overwrite order is the same on all ranks,
   * only overwrite when the rank index is greater */
  for (rank_idx = comm_rank + 1; rank_idx < blend_config->comm_size;
       rank_idx++) {
    if (blend_config->n_overlaps[rank_idx] > 0) {
      struct BlendCoords const *const rank_coords =
          blend_config->ranks_coords[rank_idx];
      AGD_TYPE const *const receive_buffer =
          blend_config->receive_buffers[rank_idx];
      int buffer_idx = 0;
      int overlap_idx;

      /* Overwrite receive buffer into blended */
      for (overlap_idx = 0; overlap_idx < blend_config->n_overlaps[rank_idx];
           overlap_idx++) {
        struct BlendCoords const *const overlap_coords =
            rank_coords + overlap_idx;
        int const channel_idx = overlap_coords->channel_idx;
        int const interval_idx = overlap_coords->interval_idx;
        size_t const interval_start = overlap_coords->interval_start;
        int const n_times = overlap_coords->n_times;
        memcpy(blended[channel_idx][interval_idx] + interval_start,
               receive_buffer + buffer_idx, (size_t)n_times * sizeof(AGD_TYPE));
        buffer_idx += n_times;
      }
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in blend_overwrite_forward_mpi\n");
  return 1;
}
#endif /* AGD_MPI */

#ifdef AGD_MPI
static int blend_mean_forward_mpi(struct BlendConfig const *const blend_config,
                                  AGD_TYPE *const *const *const blended,
                                  AGD_TYPE *const *const *const blended_count) {
  int rank_idx;

  /* Sum blended */
  if (blend_swap(blend_config, (AGD_TYPE const *const *const *)blended))
    goto err;

  for (rank_idx = 0; rank_idx < blend_config->comm_size; rank_idx++) {
    if (blend_config->n_overlaps[rank_idx] > 0) {
      struct BlendCoords const *const rank_coords =
          blend_config->ranks_coords[rank_idx];
      AGD_TYPE const *const receive_buffer =
          blend_config->receive_buffers[rank_idx];
      int buffer_idx = 0;
      int overlap_idx;

      for (overlap_idx = 0; overlap_idx < blend_config->n_overlaps[rank_idx];
           overlap_idx++) {
        struct BlendCoords const *const overlap_coords =
            rank_coords + overlap_idx;
        int const channel_idx = overlap_coords->channel_idx;
        int const interval_idx = overlap_coords->interval_idx;
        size_t const interval_start = overlap_coords->interval_start;
        int const n_times = overlap_coords->n_times;
        int time_idx;
        AGD_TYPE *const blended_interval_start =
            blended[channel_idx][interval_idx] + interval_start;
        AGD_TYPE const *const buffer_start = receive_buffer + buffer_idx;
        for (time_idx = 0; time_idx < n_times; time_idx++) {
          blended_interval_start[time_idx] += buffer_start[time_idx];
        }
        buffer_idx += n_times;
      }
    }
  }

  /* Sum count/weight */
  if (blend_swap(blend_config, (AGD_TYPE const *const *const *)blended_count))
    goto err;

  for (rank_idx = 0; rank_idx < blend_config->comm_size; rank_idx++) {
    if (blend_config->n_overlaps[rank_idx] > 0) {
      struct BlendCoords const *const rank_coords =
          blend_config->ranks_coords[rank_idx];
      AGD_TYPE const *const receive_buffer =
          blend_config->receive_buffers[rank_idx];
      int buffer_idx = 0;
      int overlap_idx;

      for (overlap_idx = 0; overlap_idx < blend_config->n_overlaps[rank_idx];
           overlap_idx++) {
        struct BlendCoords const *const overlap_coords =
            rank_coords + overlap_idx;
        int const channel_idx = overlap_coords->channel_idx;
        int const interval_idx = overlap_coords->interval_idx;
        size_t const interval_start = overlap_coords->interval_start;
        int const n_times = overlap_coords->n_times;
        int time_idx;
        AGD_TYPE *const blended_interval_start =
            blended_count[channel_idx][interval_idx] + interval_start;
        AGD_TYPE const *const buffer_start = receive_buffer + buffer_idx;
        for (time_idx = 0; time_idx < n_times; time_idx++) {
          blended_interval_start[time_idx] += buffer_start[time_idx];
        }
        buffer_idx += n_times;
      }
    }
  }

  return 0;
err:
  fprintf(stderr, "ERROR in blend_mean_forward_mpi\n");
  return 1;
}
#endif /* AGD_MPI */

static int blend_mean_forward(int const n_patches,
                              AGD_TYPE const *const *const data,
                              struct BlendConfig const *const blend_config,
                              struct BlendParams const *const blend_params,
                              int const taper_length,
                              AGD_TYPE *const *const *const blended) {
  AGD_TYPE ***blended_count = NULL;
  AGD_TYPE *taper = NULL;
  int patch_idx;
  int channel_idx;

  if (allocate_blended(blend_config, &blended_count)) goto err;
  zero_blended(blend_config, (AGD_TYPE *const *const *)blended_count);

  if (taper_length > 1) {
    int taper_idx;
    taper = (AGD_TYPE *)malloc((size_t)taper_length * sizeof(AGD_TYPE));
    if (taper == NULL) goto err;

    for (taper_idx = 0; taper_idx < taper_length; taper_idx++) {
      taper[taper_idx] = AGD_COS(AGD_HALF * (AGD_TYPE)M_PI *
                                 (AGD_TYPE)(taper_length - taper_idx) /
                                 (AGD_TYPE)(taper_length + 1));
    }
  }

  /* Sum local patches */
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    struct BlendParams const *const patch_blend_params =
        blend_params + patch_idx;
    int const trace_length = patch_blend_params->trace_length;
    int trace_idx;
    for (trace_idx = 0; trace_idx < patch_blend_params->n_traces; trace_idx++) {
      struct BlendCoords const *const trace_blend_coords =
          patch_blend_params->blend_coords + trace_idx;
      channel_idx = trace_blend_coords->channel_idx;
      if (channel_idx >= 0) {
        int const interval_idx = trace_blend_coords->interval_idx;
        size_t const interval_start = trace_blend_coords->interval_start;
        int const trace_start = trace_blend_coords->trace_start;
        int const n_times = trace_blend_coords->n_times;
        int time_idx;
        AGD_TYPE *const blended_interval_start =
            blended[channel_idx][interval_idx] + interval_start;
        AGD_TYPE *const blended_count_interval_start =
            blended_count[channel_idx][interval_idx] + interval_start;
        AGD_TYPE const *const data_start =
            data[patch_idx] + (size_t)trace_idx * (size_t)trace_length +
            (size_t)trace_start;
        for (time_idx = 0; time_idx < n_times; time_idx++) {
          int const time_to_end = n_times - 1 - time_idx;
          int const edge_dist = time_idx < time_to_end ? time_idx : time_to_end;
          int const within_taper = (edge_dist < taper_length);
          AGD_TYPE weight;
          if ((taper_length > 1) && within_taper) {
            weight = taper[edge_dist];
          } else {
            weight = AGD_ONE;
          }
          blended_interval_start[time_idx] += data_start[time_idx] * weight;
          blended_count_interval_start[time_idx] += weight;
        }
      }
    }
  }

#ifdef AGD_MPI
  /* Sum remote overlaps */
  if (blend_mean_forward_mpi(blend_config, blended,
                             (AGD_TYPE *const *const *)blended_count))
    goto err;
#endif /* AGD_MPI */

  /* Divide by weight to get mean */
  for (channel_idx = 0; channel_idx < blend_config->n_channels; channel_idx++) {
    struct ChannelIntervals const *const channel_intervals =
        blend_config->channels_intervals + channel_idx;
    int interval_idx;
    for (interval_idx = 0; interval_idx < channel_intervals->n_intervals;
         interval_idx++) {
      struct Interval const *const interval =
          channel_intervals->intervals + interval_idx;
      AGD_TYPE *const blended_interval = blended[channel_idx][interval_idx];
      AGD_TYPE const *const blended_count_interval =
          blended_count[channel_idx][interval_idx];
      size_t const n_times = (size_t)(interval->stop - interval->start);
      size_t time_idx;
      for (time_idx = 0; time_idx < n_times; time_idx++) {
        if (blended_count_interval[time_idx] > AGD_ZERO) {
          blended_interval[time_idx] /= blended_count_interval[time_idx];
        }
      }
    }
  }

  free_blended(blend_config, &blended_count);
  free(taper);
  taper = NULL;
  return 0;
err:
  fprintf(stderr, "ERROR in blend_mean_forward\n");
  free_blended(blend_config, &blended_count);
  free(taper);
  taper = NULL;
  return 1;
}

static void blend_adjoint(AGD_TYPE const *const *const *const blended,
                          struct BlendParams const *const blend_params,
                          AGD_TYPE *const data) {
  int const trace_length = blend_params->trace_length;
  int trace_idx;
  for (trace_idx = 0; trace_idx < blend_params->n_traces; trace_idx++) {
    struct BlendCoords const *const trace_blend_coords =
        blend_params->blend_coords + trace_idx;
    int const channel_idx = trace_blend_coords->channel_idx;
    if (channel_idx < 0) {
      memset(data + (size_t)trace_idx * (size_t)trace_length, 0,
             (size_t)trace_length * sizeof(AGD_TYPE));
    } else {
      int const interval_idx = trace_blend_coords->interval_idx;
      size_t const interval_start = trace_blend_coords->interval_start;
      int const trace_start = trace_blend_coords->trace_start;
      int const n_times = trace_blend_coords->n_times;
      AGD_TYPE const *const blended_interval_start =
          blended[channel_idx][interval_idx] + interval_start;
      AGD_TYPE *const data_start =
          data + (size_t)trace_idx * (size_t)trace_length + (size_t)trace_start;
      memset(data + (size_t)trace_idx * (size_t)trace_length, 0,
             (size_t)trace_start * sizeof(AGD_TYPE));
      memcpy(data_start, blended_interval_start,
             (size_t)n_times * sizeof(AGD_TYPE));
      memset(
          data_start + n_times, 0,
          (size_t)(trace_length - (trace_start + n_times)) * sizeof(AGD_TYPE));
    }
  }
}

static int blend_input(AGD_TYPE const *const *const data, int const n_patches,
                       int const *const volumes, int const *const n_dims,
                       int const *const *const data_shapes,
                       long int const *const *const shottimes,
                       int const *const *const channels,
                       enum AGDTraceType const *const *const trace_types,
                       struct BlendConfig const *const blend_config,
                       AGD_TYPE ****const blended) {
  int patch_idx;
  struct BlendParams *blend_params = (struct BlendParams *)calloc(
      (size_t)n_patches, sizeof(struct BlendParams));
  if (blend_params == NULL) goto err;

  /* Allocate blended */
  if (allocate_blended(blend_config, blended)) goto err;
  zero_blended(blend_config, (AGD_TYPE *const *const *)*blended);

  /* Set blend_params */
  for (patch_idx = 0; patch_idx < n_patches; patch_idx++) {
    int volume_idx = volumes[patch_idx];
    int const n_traces =
        get_n_traces(n_dims[volume_idx], data_shapes[patch_idx]);
    int const n_times = data_shapes[patch_idx][n_dims[volume_idx] - 1];
    if (set_blend_params(n_traces, n_times, shottimes[patch_idx],
                         channels[patch_idx], trace_types[patch_idx],
                         blend_config, blend_params + patch_idx))
      goto err;
  }

  /* Blend */
  blend_mean_forward(n_patches, data, blend_config, blend_params, 0,
                     (AGD_TYPE *const *const *)*blended);
  apply_mute(blend_config, (AGD_TYPE *const *const *)*blended);

  free_blend_params_array(n_patches, &blend_params);

  return 0;
err:
  fprintf(stderr, "ERROR in blend_input\n");
  free_blend_params_array(n_patches, &blend_params);
  return 1;
}
