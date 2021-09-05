#include <stdio.h>
#include <stdlib.h>

#ifndef EPS1
#define EPS1 0.00001f
#endif

#ifndef EPS2
#define EPS2 0.00001f
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
  size_t i;

  if (d_test == NULL || d_true == NULL) {
    fprintf(stderr, "ERROR allocating d_test or d_true\n");
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

  for (i = 0; i < N; i++) {
    if ((d_test[i] - d_true[i]) / (d_true[i] * d_true[i] + EPS1) * d_true[i] >
        EPS2) {
      fprintf(stderr, "ERROR: At element %lu, d_test = %f, d_true = %f\n", i,
              (double)d_test[i], (double)d_true[i]);
      return 1;
    }
  }

  printf("PASS\n");

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
