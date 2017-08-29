/*
 * Convolution code for 3Daudio spatializer.
 *
 * Author: Bill Gardner
 * Based on the original version by Keith Martin.
 * Copyright (c) 1995 MIT Media Lab All Rights Reserved.
 */
#include	<stdio.h>
#include	<math.h>
#include	"config.h"
#include	"sgi_fft.h"
#include	"hrtf.h"
#include	"procbuf.h"
#include	"libstuff.h"
#include	"arraylib.h"
#include	"misc.h"

#define PLAY_INPUT	FALSE	/* TRUE if just play input buffer */

/*
 * These are allocated buffers.
 */
float *in0;		/* stereo interleaved input samples, allocated */
float *in1;		/* ditto, double buffer of above, allocated */
float *out0;	/* stereo interleaved output samples, allocated */
float *out1;	/* ditto, used only for crossfades, allocated */

complex *fwd_buf;	/* buffer for forward FFT, allocated */
complex *inv_buf0;	/* buffer for inverse FFT, for left chl, allocated */
complex *inv_buf1;	/* ditto, for right chl, allocated */

float *ramp_up;		/* crossfade ramp up function (0...1) */
float *ramp_down;	/* crossfade ramp down function (1 - ramp_up) */

extern float *workspace;	/* for SGI FFT library function */

/*
 * These point to other buffers.
 */
float *p_newin;	/* ptr to new input buffer (in0 or in1) */
float *p_oldin;	/* ptr to old input buffer (in0 or in1) */

complex *filt_L;		/* left channel filter spectrum */
complex *filt_R;		/* right channel filter spectrum */
complex *old_filt_L;	/* previous left filter, when crossfading */
complex *old_filt_R;	/* previous right filter, when crossfading */

/*
 * Current gain and previous gain.
 */
double gain;
double old_gain;

void do_kdm_init(float **ibuf)
{
long i;
	in0 = f_alloc(FILT_LEN * 2);
	in1 = f_alloc(FILT_LEN * 2);
	p_newin = in0;
	p_oldin = in1;
	out0 = f_alloc(FILT_LEN * 2);
	out1 = f_alloc(FILT_LEN * 2);
	fwd_buf = c_alloc(BUF_LEN);
	inv_buf0 = c_alloc(BUF_LEN);
	inv_buf1 = c_alloc(BUF_LEN);
	get_cur_hrtf(&filt_L, &filt_R, &gain);
	/*
	 * Linear ramp table and FFT initialization.
	 */
	ramp_up = f_alloc(FILT_LEN);
	ramp_down = f_alloc(FILT_LEN);
	for (i = 0; i < FILT_LEN; i++) {
		ramp_up[i] = i / (double) FILT_LEN;
		ramp_down[i] = 1.0 - ramp_up[i];
	}
	/*
	 * First block of samples goes into p_newin.
	 */
	*ibuf = p_newin;
}

#if PLAY_INPUT

float *do_buf(float **ibuf, int flag)
{
	return *ibuf;
}

#else

/*
 * Function (should be a macro) to do spectral product, inverse
 * transform, and interleave output samples.
 */
void inv_transform(complex *filt_L, complex *filt_R, double gain, float *out)
{
float *p0, *p1, *p2;
long n;
	/*
	 * Spectral multiply and inverse transform.
	 */
	c_opt_mul(fwd_buf,filt_L,inv_buf0,BUF_LEN);
	sfft1du(1, BUF_LEN, (float *) inv_buf0, 1, workspace);
	sscal1d(FILT_LEN, gain / BUF_LEN, (float *) inv_buf0, 1);

	c_opt_mul(fwd_buf,filt_R,inv_buf1,BUF_LEN);
	sfft1du(1, BUF_LEN, (float *) inv_buf1, 1, workspace);
	sscal1d(FILT_LEN, gain / BUF_LEN, (float *) inv_buf1, 1);
	/*
	 * Interleave left and right channels into output buffer.
	 */
	p0 = (float *) inv_buf0;
	p1 = (float *) inv_buf1;
	p2 = out;
	n = FILT_LEN;
	while (n-- > 0) {
		*p2++ = *p0++;
		*p2++ = *p1++;
	}
}

float *do_buf(float **p_ibuf, int flag)
{
float *p0, *p1, *p2;
float *r0, *r1;
long n;

	if (flag) {
		old_filt_L = filt_L;
		old_filt_R = filt_R; 
		old_gain = gain;
		get_cur_hrtf(&filt_L, &filt_R, &gain);
	}

	/*
	 * Read left channel input samples into fwd_buf. Left input
	 * samples must be de-interleaved from stereo input. First
	 * read old input, then new input.
	 */
	p0 = (float *) fwd_buf;
	p1 = p_oldin;
	n = FILT_LEN;
	while (n-- > 0) {
		*p0++ = *p1;
		p1 += 2;
	}
	p1 = p_newin;
	n = FILT_LEN;
	while (n-- > 0) {
		*p0++ = *p1;
		p1 += 2;
	}
	/*
	 * Do forward FFT (real packed data).
	 */
	sfft1du(-1, BUF_LEN, (float *) fwd_buf, 1, workspace);

	/*
	 * Spectral multiply, inverse transform, and interleave
	 * output samples.
	 */
	inv_transform(filt_L, filt_R, gain, out0);

	/*
	 * If crossfading, calculate convolution with old filter
	 * and crossfade two results.
	 */
	if (flag) {
		inv_transform(old_filt_L, old_filt_R, old_gain, out1);
		/*
		 * Crossfade using pointers.
		 */
		r0 = ramp_up;
		r1 = ramp_down;
		p0 = out0;
		p1 = out1;
		p2 = out0;
		n = FILT_LEN;
		while (n-- > 0) {
			*p2++ = *p0++ * *r0 + *p1++ * *r1;
			*p2++ = *p0++ * *r0++ + *p1++ * *r1++;
		}
	}

	/*
	 * Swap input buffer pointers.
	 */
	p0 = p_newin;
	p_newin = p_oldin;
	p_oldin = p0;

	/*
	 * Change caller's pointer to input buffer.
	 */
	*p_ibuf = p_newin;

	/*
	 * Return pointer to output buffer.
	 */
	return out0;
}

#endif
