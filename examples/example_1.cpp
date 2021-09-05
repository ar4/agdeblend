#include <fstream>
#include <vector>
extern "C" {
#include "agdeblend.h"
}

int main(void) {
  int const n_patches = 1;
  int const n_traces[1] = {128};
  int const n_times[1] = {512};
  std::vector<std::vector<int long> > shottimes;
  std::vector<int long *> shottimes_ptr;
  std::vector<std::vector<int> > channels;
  std::vector<int *> channels_ptr;
  std::vector<std::vector<enum AGDTraceType> > trace_types;
  std::vector<enum AGDTraceType *> trace_types_ptr;
  std::vector<std::vector<float> > data;
  std::vector<float *> data_ptr;
  enum AGDBlendMode blend_mode = AGDBlendSum;
  int taper_length = 0;
  std::ifstream infile;
  std::ofstream outfile;

  /* Allocate memory */
  shottimes.resize(n_patches);
  shottimes[0].resize(n_traces[0]);
  channels.resize(n_patches);
  channels[0].resize(n_traces[0]);
  trace_types.resize(n_patches);
  trace_types[0].resize(n_traces[0]);
  data.resize(n_patches);
  data[0].resize(n_traces[0] * n_times[0]);

  /* Set pointers */
  shottimes_ptr.push_back(&shottimes[0][0]);
  channels_ptr.push_back(&channels[0][0]);
  trace_types_ptr.push_back(&trace_types[0][0]);
  data_ptr.push_back(&data[0][0]);

  /* Load from input files */
  infile.open("data_1_shottimes.bin",
              std::ifstream::in | std::ifstream::binary);
  infile.read(reinterpret_cast<char *>(shottimes[0].data()),
              sizeof(long int) * n_traces[0]);
  infile.close();

  infile.open("data_1_channels.bin", std::ifstream::in | std::ifstream::binary);
  infile.read(reinterpret_cast<char *>(channels[0].data()),
              sizeof(int) * n_traces[0]);
  infile.close();

  infile.open("data_1_trace_types.bin",
              std::ifstream::in | std::ifstream::binary);
  infile.read(reinterpret_cast<char *>(trace_types[0].data()),
              sizeof(enum AGDTraceType) * n_traces[0]);
  infile.close();

  infile.open("data_1_true_data.bin",
              std::ifstream::in | std::ifstream::binary);
  infile.read(reinterpret_cast<char *>(data[0].data()),
              sizeof(int) * n_traces[0] * n_times[0]);
  infile.close();

  /* Blend */
  agd_blend(n_patches, n_traces, n_times, shottimes_ptr.data(),
            channels_ptr.data(), trace_types_ptr.data(), data_ptr.data(),
            blend_mode, taper_length, n_patches, n_traces, n_times,
            shottimes_ptr.data(), channels_ptr.data(), trace_types_ptr.data(),
            data_ptr.data());

  /* Write to file */
  outfile.open("out/data_1_blended_data.bin",
               std::ofstream::out | std::ofstream::binary);
  outfile.write(reinterpret_cast<char *>(data[0].data()),
                sizeof(int) * n_traces[0] * n_times[0]);
  outfile.close();

  return 0;
}
