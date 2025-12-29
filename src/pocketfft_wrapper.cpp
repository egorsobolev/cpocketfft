extern "C"
{
#include "fft.h"
}

#include "pocketfft/pocketfft_hdronly.h"
#include <stdio.h>

using std::size_t;
using pocketfft::shape_t;
using pocketfft::stride_t;


void
fft_make_shape (size_t n, size_t *nn, size_t itemsize_real, size_t itemsize_cmplx,
                shape_t & shape, stride_t & stride_real, stride_t & stride_cmplx,
                shape_t & axes)
{
	size_t k, i, j;

	i = n - 1;
	k = shape[i] = nn[i];
	stride_real[i] = itemsize_real;
	itemsize_real *= k;
	stride_cmplx[i] = itemsize_cmplx;
	itemsize_cmplx *= (k / 2) + 1;
	axes[i] = i;
	for (j = 1; j < n; j++) {
		i = n - 1 - j;
		k = shape[i] = nn[i];
		stride_real[i] = itemsize_real;
		itemsize_real *= k;
		stride_cmplx[i] = itemsize_cmplx;
		itemsize_cmplx *= k;
		axes[i] = i;
	}
}


void
fft_r2c (size_t n, size_t *nn, const double *in, complex_double * out)
{
	int i;
	shape_t shape (n), axes (n);
	stride_t stride_in (n), stride_out (n);
	fft_make_shape (n, nn, sizeof (double), sizeof (complex_double),
	                shape, stride_in, stride_out, axes);
	pocketfft::r2c (shape, stride_in, stride_out, axes, pocketfft::FORWARD,
	                in, reinterpret_cast < std::complex <double >*>(out),
	                1.0, 0);
}

void
fftf_r2c (size_t n, size_t *nn, const float *in, complex_float * out)
{
	shape_t shape (n), axes (n);
	stride_t stride_in (n), stride_out (n);
	fft_make_shape (n, nn, sizeof (float), sizeof (complex_float),
	                shape, stride_in, stride_out, axes);
	pocketfft::r2c (shape, stride_in, stride_out, axes, pocketfft::FORWARD,
	                in, reinterpret_cast < std::complex <float >*>(out),
	                1.0f, 0);
}

void
fft_c2r (size_t n, size_t *nn, const complex_double * in, double *out)
{
	int i;
	shape_t shape (n), axes (n);
	stride_t stride_in (n), stride_out (n);
	fft_make_shape (n, nn, sizeof (double), sizeof (complex_double),
	                shape, stride_out, stride_in, axes);
	pocketfft::c2r (shape, stride_in, stride_out, axes, pocketfft::BACKWARD,
	                reinterpret_cast < const std::complex <double >*>(in), out,
	                1.0, 0);
}

void
fftf_c2r (size_t n, size_t *nn, const complex_float * in, float *out)
{
	shape_t shape (n), axes (n);
	stride_t stride_in (n), stride_out (n);
	fft_make_shape (n, nn, sizeof (float), sizeof (complex_float),
	                shape, stride_out, stride_in, axes);
	pocketfft::c2r (shape, stride_in, stride_out, axes, pocketfft::BACKWARD,
	                reinterpret_cast < const std::complex <float >*>(in), out,
	                1.0f, 0);
}
