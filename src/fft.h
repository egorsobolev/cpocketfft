#ifndef CPOCKETFFT_FFT_H
#define CPOCKETFFT_FFT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct CMPLX_DBL
{
	double r, i;
} complex_double;

typedef struct CMPLX_FLT
{
	float r, i;
} complex_float;


void fft_r2c(size_t n, size_t *shape, const double *in, complex_double *out,
             double scale);
void fftf_r2c(size_t n, size_t *shape, const float *in, complex_float *out,
              float scale);
void fft_c2r(size_t n, size_t *shape, const complex_double *in, double *out,
             double scale);
void fftf_c2r(size_t n, size_t *shape, const complex_float *in, float *out,
              float scale);

#ifdef __cplusplus
}
#endif

#endif
