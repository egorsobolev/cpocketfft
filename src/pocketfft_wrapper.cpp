extern "C"
{
#include "fft.h"
}

#include "pocketfft/pocketfft_hdronly.h"
#include <stdlib.h>


using std::size_t;
using pocketfft::shape_t;
using pocketfft::stride_t;


void
fft_make_shape_r2c (size_t n, size_t *nn, size_t itemsize_real, size_t itemsize_cmplx,
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
fft_make_shape (size_t n, size_t *nn, size_t itemsize, shape_t & shape,
                stride_t & stride, shape_t & axes)
{
	size_t k, i, j;

	for (j = 0; j < n; j++) {
		i = n - 1 - j;
		k = shape[i] = nn[i];
		stride[i] = itemsize;
		itemsize *= k;
		axes[i] = i;
	}
}

int cpocketfft_num_threads = 0;

void
cpocketfft_set_num_threads(int num_threads)
{
	cpocketfft_num_threads = num_threads;
}

void
fft_r2c (size_t n, size_t *nn, const double *in, complex_double * out, double scale)
{
	int i;
	shape_t shape (n), axes (n);
	stride_t stride_in (n), stride_out (n);
	fft_make_shape_r2c (n, nn, sizeof (double), sizeof (complex_double),
	                    shape, stride_in, stride_out, axes);
	pocketfft::r2c (shape, stride_in, stride_out, axes, pocketfft::FORWARD,
	                in, reinterpret_cast < std::complex <double >*>(out),
	                scale, cpocketfft_num_threads);
}

void
fftf_r2c (size_t n, size_t *nn, const float *in, complex_float * out, float scale)
{
	shape_t shape (n), axes (n);
	stride_t stride_in (n), stride_out (n);
	fft_make_shape_r2c (n, nn, sizeof (float), sizeof (complex_float),
                        shape, stride_in, stride_out, axes);
	pocketfft::r2c (shape, stride_in, stride_out, axes, pocketfft::FORWARD,
	                in, reinterpret_cast < std::complex <float >*>(out),
	                scale, cpocketfft_num_threads);
}

void
fft_c2r (size_t n, size_t *nn, const complex_double * in, double *out,
         double scale)
{
	int i;
	shape_t shape (n), axes (n);
	stride_t stride_in (n), stride_out (n);
	fft_make_shape_r2c (n, nn, sizeof (double), sizeof (complex_double),
	                    shape, stride_out, stride_in, axes);
	pocketfft::c2r (shape, stride_in, stride_out, axes, pocketfft::BACKWARD,
	                reinterpret_cast < const std::complex <double >*>(in), out,
	                scale, cpocketfft_num_threads);
}

void
fftf_c2r (size_t n, size_t *nn, const complex_float * in, float *out, float scale)
{
	shape_t shape (n), axes (n);
	stride_t stride_in (n), stride_out (n);
	fft_make_shape_r2c (n, nn, sizeof (float), sizeof (complex_float),
	                    shape, stride_out, stride_in, axes);
	pocketfft::c2r (shape, stride_in, stride_out, axes, pocketfft::BACKWARD,
	                reinterpret_cast < const std::complex <float >*>(in), out,
	                scale, cpocketfft_num_threads);
}

void
fft_dst(size_t n, size_t *nn, int type, const double *in, double *out,
        double scale, int ortho)
{
	shape_t shape (n), axes (n);
	stride_t stride (n);
	fft_make_shape(n, nn, sizeof(double), shape, stride, axes);
	pocketfft::dst(shape, stride, stride, axes, type, in, out,
	               scale, ortho, cpocketfft_num_threads);
}

void
fftf_dst(size_t n, size_t *nn, int type, const float *in, float *out,
         float scale, int ortho)
{
	shape_t shape (n), axes (n);
	stride_t stride (n);
	fft_make_shape(n, nn, sizeof(float), shape, stride, axes);
	pocketfft::dst(shape, stride, stride, axes, type, in, out,
	               scale, ortho, cpocketfft_num_threads);
}

void
fft_dct(size_t n, size_t *nn, int type, const double *in, double *out,
        double scale, int ortho)
{
	shape_t shape (n), axes (n);
	stride_t stride (n);
	fft_make_shape(n, nn, sizeof(double), shape, stride, axes);
	pocketfft::dct(shape, stride, stride, axes, type, in, out,
	               scale, ortho, cpocketfft_num_threads);
}

void
fftf_dct(size_t n, size_t *nn, int type, const float *in, float *out,
         float scale, int ortho)
{
	shape_t shape (n), axes (n);
	stride_t stride (n);
	fft_make_shape(n, nn, sizeof(float), shape, stride, axes);
	pocketfft::dct(shape, stride, stride, axes, type, in, out,
	               scale, ortho, cpocketfft_num_threads);
}
