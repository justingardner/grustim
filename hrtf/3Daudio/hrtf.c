/*
 * Code to load/transform/lookup HRTFs and maintain MIDI control
 * values for 3Daudio spatializer.
 *
 * Author: Bill Gardner
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

/*************************************************************
 * Default root path, can be overridden by setting environment
 * variable HRTFROOT.
 */
#define ROOT_ENV	"HRTFROOT"
#define ROOT_PATH	"/ti/u/billg/hrtf/diffuse32k"
#define PATH_FMT	"%s/elev%d/H%de%03da%s"
#define EXTENSION	".res"
#define HRTF_LEN	256

char rootpath[256];

/*
 * An HRTF is stored as left and right spectra for easy convolution.
 * The compact HRTFs are 128 time samples, which give us 128 complex
 * spectral points (which are conjugate symmetric). Depending on the
 * type of fft used, we will either have 128 complex points (fft.h)
 * or 65 complex points (sgi_fft.h).
 */
typedef struct {
	complex *left;
	complex *right;
} HRTF_DATUM;

/*
 * A two-dimensional array of HRTF_DATA, stored by elevation
 * and azimuth. Use get_indices() to get the indices for this
 * array.
 */
HRTF_DATUM **hrtf_data;

#define MIN_ELEV	-40
#define MAX_ELEV	90
#define ELEV_INC	10
#define N_ELEV		(((MAX_ELEV - MIN_ELEV) / ELEV_INC) + 1)

#define MIN_AZIM	-180
#define MAX_AZIM	180

#define MIN_ATTEN	0	/* in dB */
#define MAX_ATTEN	20

/*
 * This array gives the total number of azimuths measured
 * per elevation, and hence the AZIMUTH INCREMENT. Note that
 * only azimuths up to and including 180 degrees actually
 * exist in file system (because the data is symmetrical.
 */
int elev_data[N_ELEV] = {
56, 60, 72, 72, 72, 72, 72, 60, 56, 45, 36, 24, 12, 1 };

/*
 * MIDI control stuff.
 */
#define ELEV_CTL	1
#define AZIM_CTL	2
#define ATTEN_CTL	3

#define N_CTLS		64
#define MAX_CTL_VAL	127

#define CTL_VAL(ctl_num)	((double) midi_ctls[ctl_num] / MAX_CTL_VAL)

/*
 * Array of MIDI controls.
 */
int midi_ctls[N_CTLS];

/*
 * Workspace for SGI FFT library. Used by hrtf.c and procbuf.c.
 */
float *workspace;


/**************************************************************************
 * Functions
 */


static int round(double val)
{
	return (val > 0) ? val + 0.5 : val - 0.5;
}

/*
 * Return the number of azimuths actually stored in file system.
 */
int get_nfaz(int el_index)
{
	return (elev_data[el_index] / 2) + 1;
}

/*
 * Get (closest) elevation index for given elevation.
 * Elevation is clipped to legal range.
 */
int get_el_index(double elev)
{
int el_index;
	el_index = round((elev - MIN_ELEV) / ELEV_INC);
	if (el_index < 0) el_index = 0;
	else if (el_index >= N_ELEV) el_index = N_ELEV - 1;
	return el_index;
}

/*
 * For a given elevation and azimuth in degrees, return the
 * indices for the proper HRTF. *p_flip will be set TRUE if
 * left and right channels need to be flipped.
 */
void get_indices(double elev, double azim, int *p_el, int *p_az, int *p_flip)
{
int naz, nfaz;
int el_index;
int az_index;
	el_index = get_el_index(elev);
	naz = elev_data[el_index];
	nfaz = get_nfaz(el_index);
	/*
	 * Coerce azimuth into legal range and calculate
	 * flip if any.
	 */
	azim = fmod(azim, 360.0);
	if (azim < 0) azim += 360;
	if (azim > 180) {
		azim = 360 - azim;
		*p_flip = TRUE;
	}
	else *p_flip = FALSE;

	/*
	 * Now 0 <= azim <= 180. Calculate index and clip to
	 * legal range just to be sure.
	 */
	az_index = round(azim / (360.0 / naz));
	if (az_index < 0) az_index = 0;
	else if (az_index >= nfaz) az_index = nfaz - 1;

	*p_el = el_index;
	*p_az = az_index;
}

/*
 * Convert indices to angles.
 */
int index_to_elev(int el_index)
{
	return MIN_ELEV + el_index * ELEV_INC;
}

int index_to_azim(int el_index, int az_index)
{
	return round(az_index * 360.0 / elev_data[el_index]);
}

/*
 * Return pathname of HRTF specified by indices.
 */
char *hrtf_name(int el_index, int az_index)
{
static char buf[128];
int elev;
int azim;

	elev = index_to_elev(el_index);
	azim = index_to_azim(el_index, az_index);

	sprintf(buf,PATH_FMT,
		rootpath,elev,elev,azim,EXTENSION);
	return buf;
}

/*
 * Compute spectrum of short array as needed by convolution alg.
 * Filter response is copied into FILT_LEN size buffer, PRE-padded
 * with FILT_LEN zeros, and transformed. The resulting size of
 * returned spectrum depends on which fft package we're using.
 */
complex *xform_chan(short *sx, long n)
{
long i;
float *x;

	x = (float *) c_alloc(BUF_LEN / 2 + 1);
	for (i = 0; i < n; i++)
		x[FILT_LEN + i] = (double) sx[i] / 32768;
	sfft1du(-1, BUF_LEN, x, 1, workspace);

	return (complex *) x;
}

/*
 * Read and transform the HRTF specified by the indices.
 * Store into pre-allocated hd. We always store into
 * 128-pt result, even if HRTF data is smaller.
 */
void read_hrtf(int el_index, int az_index, HRTF_DATUM *hd)
{
char *name;
short *buf;
short *sl, *sr;
long n, nr;
	name = hrtf_name(el_index, az_index);
	buf = file_load(name,&nr);
	n = nr / 2;
	if (n > FILT_LEN) {
		fprintf(stderr,"unexpected length of HRTF, %ld sample frames\n",n);
		exit(0);
	}
	sl = s_alloc(n);
	sr = s_alloc(n);
	s_split(buf, n, sl, sr);
	free(buf);
	/*
	 * Transform channels.
	 */
	hd->left = xform_chan(sl, n);
	hd->right = xform_chan(sr, n);
	free(sl);
	free(sr);
}

/*
 * Read and transform all the HRTFs.
 */
void read_hrtfs()
{
int el_index;
int az_index;
int nfaz;
char *root;
	root = getenv(ROOT_ENV);
	strcpy(rootpath, root ? root : ROOT_PATH);
	printf("Loading HRTFs from '%s'\n",rootpath);
	hrtf_data = (HRTF_DATUM **) ralloc(N_ELEV * sizeof(HRTF_DATUM *));
	for (el_index = 0; el_index < N_ELEV; el_index++) {
		printf("elevation %d...\n", index_to_elev(el_index));
		nfaz = get_nfaz(el_index);
		hrtf_data[el_index] = 
			(HRTF_DATUM *) ralloc(nfaz * sizeof(HRTF_DATUM));
		for (az_index = 0; az_index < nfaz; az_index++) {
			read_hrtf(el_index, az_index,
				&hrtf_data[el_index][az_index]);
		}
	}
}

void free_hrtfs()
{
int el_index, az_index;
int nfaz;
	for (el_index = 0; el_index < N_ELEV; el_index++) {
		nfaz = get_nfaz(el_index);
		for (az_index = 0; az_index < nfaz; az_index++) {
			free(hrtf_data[el_index][az_index].left);
			free(hrtf_data[el_index][az_index].right);
		}
		free(hrtf_data[el_index]);
	}
	free(hrtf_data);
}

/*
 * Called with MIDI control numbers and values.
 */
void do_ctl(int ctl_num, int ctl_val)
{
	if (ctl_num >= 0 && ctl_num < N_CTLS)
		midi_ctls[ctl_num] = ctl_val;
}

#if TESTING
/*
 * Code for testing functionality. Set up stub routines for
 * get_hrtf() ...etc, that just return simple filter responses.
 */
complex *test_left;
complex *test_right;

void setup_test_hrtfs()
{
short *sx;
long nx;
	printf("setup test hrtfs...\n");
	nx = FILT_LEN;
	sx = s_alloc(nx);
	sx[2] = 32767;
	test_left = xform_chan(sx, nx);
	printf("test_left half spectrum:\n");
	c_print(test_left, FILT_LEN + 1);
	s_clear(sx,nx);
	sx[0] = 16384;
	sx[1] = -16384;
	test_right = xform_chan(sx, nx);
	printf("test_right half spectrum:\n");
	c_print(test_right, FILT_LEN + 1);
}

void get_hrtf(double elev, double azim, complex **p_left, complex **p_right)
{
	printf("get_hrtf...\n");
	*p_left = test_left;
	*p_right = test_right;
}

void get_cur_hrtf(complex **p_left, complex **p_right, double *p_gain)
{
	printf("get_cur_hrtf...\n");
	*p_left = test_left;
	*p_right = test_right;
	*p_gain = 1.0;
}

#else

/*
 * Current and last values for position. Need these to print
 * position values only when things change. Pretty gross.
 */
int cur_el_index;
int cur_az_index;
int cur_flip_flag;
int last_el_index;
int last_az_index;
int last_flip_flag;
double cur_atten;
double last_atten;
int changed;

/*
 * Get the closest HRTF to the specified elevation and azimuth
 * in degrees. Return via left and right channel pointers.
 */
void get_hrtf(double elev, double azim, complex **p_left, complex **p_right)
{
HRTF_DATUM *hd;
	/*
	 * Clip angles and convert to indices.
	 */
	get_indices(elev, azim, &cur_el_index, &cur_az_index, &cur_flip_flag);

	if (cur_el_index != last_el_index ||
		cur_az_index != last_az_index ||
		cur_flip_flag != last_flip_flag)
			changed = TRUE;
	last_el_index = cur_el_index;
	last_az_index = cur_az_index;
	last_flip_flag = cur_flip_flag;

	/*
	 * Get data and flip channels if necessary.
	 */
	hd = &hrtf_data[cur_el_index][cur_az_index];
	if (cur_flip_flag) {
		*p_left = hd->right;
		*p_right = hd->left;
	}
	else {
		*p_left = hd->left;
		*p_right = hd->right;
	}
}

/*
 * Call to get current left and right HRTFs based on
 * MIDI control values.
 */
void get_cur_hrtf(complex **p_left, complex **p_right, double *p_gain)
{
double elev, azim;
	changed = FALSE;
	elev = MIN_ELEV + (MAX_ELEV - MIN_ELEV) * CTL_VAL(ELEV_CTL);
	azim = MIN_AZIM + (MAX_AZIM - MIN_AZIM) * CTL_VAL(AZIM_CTL);
	get_hrtf(elev, azim, p_left, p_right);

	cur_atten = MIN_ATTEN + (MAX_ATTEN - MIN_ATTEN) * CTL_VAL(ATTEN_CTL);
	*p_gain = pow(10.0, -cur_atten / 20);
	if (cur_atten != last_atten) changed = TRUE;
	last_atten = cur_atten;

	if (changed) {
		printf("elev %d azim %d atten %.0lf\n",
			index_to_elev(cur_el_index),
			cur_flip_flag ? -index_to_azim(cur_el_index, cur_az_index) :
			index_to_azim(cur_el_index, cur_az_index),
			cur_atten);
	}
}

#endif

/*
 * Initialization function.
 */
void do_init(float **p_ibuf)
{
	workspace = sfft1dui(BUF_LEN, NULL);
#if TESTING
	setup_test_hrtfs();
#else
	read_hrtfs();
#endif
	do_kdm_init(p_ibuf);
}
