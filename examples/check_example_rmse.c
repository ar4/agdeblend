#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef EPS
#define EPS 0.00001f
#endif

#ifndef MAX_ERROR_FRAC
#define MAX_ERROR_FRAC (1.0f / 50.0f)
#endif

#ifndef FILENAME
#error Please define FILENAME
#endif

#ifndef OUT_FILENAME
#define OUT_FILENAME FILENAME
#endif

#ifndef N
#error Please define N
#endif

int main(void) {
  FILE *fid_test = fopen("out/" OUT_FILENAME, "rb");
  FILE *fid_true = fopen(FILENAME, "rb");
  float *d_test = (float *)malloc(N * sizeof(float));
  float *d_true = (float *)malloc(N * sizeof(float));
  float true_rms = 0.0f;
  float error_rms = 0.0f;
  size_t i;

  if (d_test == NULL || d_true == NULL) {
    fprintf(stderr, "ERROR allocating d_test or d_true\n");
    goto err;
  }

  if (fid_test == NULL) {
    fprintf(stderr, "ERROR opening %s\n", "out/" OUT_FILENAME);
    goto err;
  }

  if (fid_true == NULL) {
    fprintf(stderr, "ERROR opening %s\n", FILENAME);
    goto err;
  }

  i = fread(d_test, sizeof(float), N, fid_test);

  if (i != N) {
    fprintf(stderr, "ERROR reading d_test\n");
    goto err;
  }

  i = fread(d_true, sizeof(float), N, fid_true);

  if (i != N) {
    fprintf(stderr, "ERROR reading d_true\n");
    goto err;
  }

  /* Calculate the RMS amplitude of the true data */
  for (i = 0; i < N; i++) {
    true_rms += d_true[i] * d_true[i];
  }
  true_rms /= (float)N;
  true_rms = sqrtf(true_rms);

  /* Calculate the RMS amplitude of the error */
  for (i = 0; i < N; i++) {
    float const error = d_true[i] - d_test[i];
    error_rms += error * error;
  }
  error_rms /= (float)N;
  error_rms = sqrtf(error_rms);

  if (error_rms / (true_rms + EPS) > MAX_ERROR_FRAC) {
    fprintf(stderr, "ERROR: RMS error fraction = %f\n",
            error_rms / (true_rms + EPS));
    return 1;
  }

  printf("PASS %f\n", error_rms / (true_rms + EPS));

  fclose(fid_test);
  fclose(fid_true);
  free(d_test);
  free(d_true);
  return 0;

err:
  fclose(fid_test);
  fclose(fid_true);
  free(d_test);
  free(d_true);
  return 1;
}
