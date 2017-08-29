#include <math.h>
#include <sys/types.h>
#include <malloc.h>
#include <stdlib.h>

#ifndef _SGI_FFT_
#define	_SGI_FFT_

#define	    FACTOR_SPACE    15

/* *******************************************************
    Complex structures definitions
******************************************************* */

typedef struct {
    float re;
    float im;
} complex;

typedef struct {
    double re;
    double im;
} zomplex;


/* *******************************************************
    C Functions prototypes
******************************************************* */
/* *******************************************************
	complex <---> complex FFTs
******************************************************* */
complex *cfft1di( int n, complex *save);
int cfft1d( int job, int n, complex *array, int inc, complex *save);

complex *cfft2di( int n1, int n2, complex *save);
int cfft2d( int job, int n1, int n2, complex *array, int ld, complex *save);

complex *cfft3di( int n1, int n2, int n3, complex *save);
int cfft3d( int job, int n1, int n2, int n3, complex *array, int ld1, int ld2, complex *save);


/* *******************************************************
	zomplex <---> zomplex FFTs
******************************************************* */
zomplex *zfft1di( int n, zomplex *save);
int zfft1d( int job, int n, zomplex *array, int inc, zomplex *save);

zomplex *zfft2di( int n1, int n2, zomplex *save);
int zfft2d( int job, int n1, int n2, zomplex *array, int ld, zomplex *save);

zomplex *zfft3di( int n1, int n2, int n3, zomplex *save);
int zfft3d( int job, int n1, int n2, int n3, zomplex *array, int ld1, int ld2, zomplex *save);


/* *******************************************************
	real <---> complex FFTs
******************************************************* */
float *sfft1di( int n, float *save);
#define sfft1dui    sfft1di
int sfft1d( int job, int n, float *array, int inc, float *save);
int sfft1du( int job, int n, float *array, int inc, float *save);

float *sfft2di( int n1, int n2, float *save);
float *sfft2dui( int n1, int n2, float *save);
int sfft2d( int job, int n1, int n2, float *array, int ld, float *save);
int sfft2du( int job, int n1, int n2, float *array, int ld, float *save);

float *sfft3di( int n1, int n2, int n3, float *save);
float *sfft3dui( int n1, int n2, int n3, float *save);
int sfft3d( int job, int n1, int n2, int n3, float *array, int ld1, int ld2, float *save);
int sfft3du( int job, int n1, int n2, int n3, float *array, int ld1, int ld2, float *save);


/* *******************************************************
	double <---> zomplex FFTs
******************************************************* */
double *dfft1di( int n, double *save);
#define dfft1dui    dfft1di
int dfft1d( int job, int n, double *array, int inc, double *save);
int dfft1du( int job, int n, double *array, int inc, double *save);

double *dfft2di( int n1, int n2, double *save);
double *dfft2dui( int n1, int n2, double *save);
int dfft2d( int job, int n1, int n2, double *array, int ld, double *save);
int dfft2du( int job, int n1, int n2, double *array, int ld, double *save);

double *dfft3di( int n1, int n2, int n3, double *save);
double *dfft3dui( int n1, int n2, int n3, double *save);
int dfft3d( int job, int n1, int n2, int n3, double *array, int ld1, int ld2, double *save);
int dfft3du( int job, int n1, int n2, int n3, double *array, int ld1, int ld2, double *save);


/* *******************************************************
    Fortran Subroutines prototypes
******************************************************* */
/* *******************************************************
	complex <---> complex FFTs
******************************************************* */
void cfft1di_( int *n, complex *save);
int cfft1d_( int *job, int *n, complex *array, int *inc, complex *save);

void cfft2di_( int *n1, int *n2, complex *save);
void cfft2d_( int *job, int *n1, int *n2, complex *array, int *ld, complex *save);

void cfft3di_( int *n1, int *n2, int *n3, complex *save);
void cfft3d_( int *job, int *n1, int *n2, int *n3, complex *array, int *ld1, int *ld2, complex *save);

/* ****************************
	zomplex <---> zomplex FFTs
**************************** */
void zfft1di_( int *n, zomplex *save);
void zfft1d_( int *job, int *n, zomplex *array, int *inc, zomplex *save);

void zfft2di_( int *n1, int *n2, zomplex *save);
void zfft2d_( int *job, int *n1, int *n2, zomplex *array, int *ld, zomplex *save);

void zfft3di_( int *n1, int *n2, int *n3, zomplex *save);
void zfft3d_( int *job, int *n1, int *n2, int *n3, zomplex *array, int *ld1, int *ld2, zomplex *save);

/* ********************************************************
	real <---> complex FFTs
******************************************************** */
void *sfft1di_( int *n, float *save);
#define sfft1dui_ sfft1di_
int *sfft1d_( int *job, int *n, float *array, int *inc, float *save);
int *sfft1du_( int *job, int *n, float *array, int *inc, float *save);

void *sfft2di_( int *n1, int *n2, void *save);
void *sfft2dui_( int *n1, int *n2, void *save);
int *sfft2d_( int *job, int *n1, int *n2, void *array, int *ld, void *save);
int *sfft2du_( int *job, int *n1, int *n2, void *array, int *ld, void *save);

void *sfft3di_( int *n1, int *n2, int *n3, void *save);
void *sfft3dui_( int *n1, int *n2, int *n3, void *save);
int *sfft3d_( int *job, int *n1, int *n2, int *n3, void *array, int *ld1, int *ld2, void *save);
int *sfft3du_( int *job, int *n1, int *n2, int *n3, void *array, int *ld1, int *ld2, void *save);

/* ********************************************************
	double <---> zomplex FFTs
******************************************************** */
void *dfft1di_( int *n, double *save);
#define dfft1dui_ dfft1di_
int *dfft1d_( int *job, int *n, double *array, int *inc, double *save);
int *dfft1du_( int *job, int *n, double *array, int *inc, double *save);

void *dfft2di_( int *n1, int *n2, void *save);
void *dfft2dui_( int *n1, int *n2, void *save);
int *dfft2d_( int *job, int *n1, int *n2, void *array, int *ld, void *save);
int *dfft2du_( int *job, int *n1, int *n2, void *array, int *ld, void *save);

void *dfft3di_( int *n1, int *n2, int *n3, void *save);
void *dfft3dui_( int *n1, int *n2, int *n3, void *save);
int *dfft3d_( int *job, int *n1, int *n2, int *n3, void *array, int *ld1, int *ld2, void *save);
int *dfft3du_( int *job, int *n1, int *n2, int *n3, void *array, int *ld1, int *ld2, void *save);

/************************************************************************
    Product modules ... 
	    Performs convolution in 1 Domain by Product in the other
********************************************************************** */

void cprod1d( int n, complex *y, int incy, complex *filter, int incx);
void zprod1d( int n, complex *y, int incy, complex *filter, int incx);
void sprod1du( int n, complex *y, int incy, complex *filter, int incx);
void dprod1du( int n, complex *y, int incy, complex *filter, int incx);
void sprod1d( int n, complex *y, int incy, complex *filter, int incx);
void dprod1d( int n, complex *y, int incy, complex *filter, int incx);
void csprod1d( int n, complex *y, int iry, int incy, complex *filter, 
		int irx, int incx);
void zdprod1d( int n, complex *y, int iry, int incy, complex *filter, 
		int irx, int incx);
		
void cprod2d( int n1, int n2, complex *y, int ldy, complex *filter, int ldx);
void zprod2d( int n1, int n2, zomplex *y, int ldy, zomplex *filter, int ldx);
void sprod2du( int n1, int n2, float *y, int ldy, float *filter, int ldx);
void dprod2du( int n1, int n2, double *y, int ldy, double *filter, int ldx);
void sprod2d( int n1, int n2, float *y, int ldy, float *filter, int ldx);
void dprod2d( int n1, int n2, double *y, int ldy, double *filter, int ldx);

void cprod3d( int n1, int n2, int n3, complex *y, int ldy1, int ldy2, 
	    complex *filter, int ldx1, int ldx2);
void zprod3d( int n1, int n2, int n3, complex *y, int ldy1, int ldy2, 
	    complex *filter, int ldx1, int ldx2);
void sprod3du( int n1, int n2, int n3, complex *y, int ldy1, int ldy2, 
	    complex *filter, int ldx1, int ldx2);
void dprod3du( int n1, int n2, int n3, complex *y, int ldy1, int ldy2, 
	    complex *filter, int ldx1, int ldx2);
void sprod3d( int n1, int n2, int n3, complex *y, int ldy1, int ldy2, 
	    complex *filter, int ldx1, int ldx2);
void dprod3d( int n1, int n2, int n3, complex *y, int ldy1, int ldy2, 
	    complex *filter, int ldx1, int ldx2);
/************************************************************************
    Scaling modules ... 
	    Scale the sequence by value APLHA ... 
		to keep absolute values after Direct+Inverse transform.
********************************************************************** */

void cscal1d( int n, float alpha, complex *y, int inc);
void zscal1d( int n, double alpha, zomplex *y, int inc);
void sscal1d( int n, float alpha, float *y, int inc);
void dscal1d( int n, double alpha, double *y, int inc);

void cscal2d( int nx, int ny, float alpha, complex *y, int ld);
void zscal2d( int nx, int ny, double alpha, zomplex *y, int ld);
void sscal2d( int nx, int ny, float alpha, float *y, int ld);
void dscal2d( int nx, int ny, double alpha, double *y, int ld);

void cscal3d( int nx, int ny, int nz, float alpha, complex *y, int ld1,int ld2);
void zscal3d( int nx, int ny, int nz, double alpha, zomplex *y,int ld1,int ld2);
void sscal3d( int nx, int ny, int nz, float alpha, float *y, int ld1,int ld2);
void dscal3d( int nx, int ny, int nz, double alpha, double *y, int ld1,int ld2);

#endif
