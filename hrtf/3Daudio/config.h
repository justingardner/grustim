/*
 * Author: Bill Gardner
 * Copyright (c) 1995 MIT Media Lab All Rights Reserved.
 */

#define TESTING		0		/* TRUE if just testing convolution */

/*
 * Length of the HRTF filters (rounded up to FFT size).
 */
#if TESTING
#define FILT_LEN	8L
#else
#define FILT_LEN	128L
#endif

/*
 * FFT buffers are this many complex points.
 */
#define BUF_LEN		(FILT_LEN * 2)
