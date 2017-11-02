/*
 * Miscellaneous routines for spatializer demo taken from
 * various library modules.
 *
 * Author: Bill Gardner
 * Copyright (c) 1995 MIT Media Lab All Rights Reserved.
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	"complex.h"
#include	"arraylib.h"
#include	"libstuff.h"

void s_split(
	short *stereo,	/* ptr to stereo samples */
	long ns,		/* number stereo frames */
	short *left,	/* left result, pre-allocated */
	short *right)	/* right result, pre-allocated */
{
short *end;
	end = stereo + ns * 2;
	while (stereo < end) {
		*left++ = *stereo++;
		*right++ = *stereo++;
	}
}

static complex cMul(a,b)
complex a, b;
{
complex c;
	c.re = a.re * b.re - a.im * b.im;
	c.im = a.im * b.re + a.re * b.im;
	return c;
}

/*
 * Optimized multiply of two conjugate symmetric arrays (c = a * b).
 */
void c_opt_mul(a,b,c,n)
complex *a, *b, *c;
long n;
{
register long i, half = n >> 1;
	for (i = 0; i <= half; i++)
		c[i] = cMul(*a++,*b++);
	for (i = 1; i < half; i++) {
		c[half+i].re = c[half-i].re;
		c[half+i].im = -c[half-i].im;
	}
}

static long get_file_size(fp)
FILE *fp;
{
long n;
	fseek(fp,0L,2);
	n = ftell(fp);
	fseek(fp,0L,0);
	return n;
}

#define BLOCK_SIZE	1024
#define MIN(a,b)	(((a) < (b)) ? (a) : (b))

short *file_load(fname,p_n)
char *fname;
long *p_n;
{
FILE *fp;
long n;
int nr;
short *buf;
short *buf0;
long count;
	if ((fp = fopen(fname,"rb")) == NULL) {
		fprintf(stderr,"can't fopen %s\n",fname);
		exit(0);
	}
	*p_n = n = get_file_size(fp) / sizeof(short);
	buf0 = buf = s_alloc(n);
	count = n;
	while (count > 0) {
		nr = MIN(BLOCK_SIZE,count);
		if (fread(buf,sizeof(short),nr,fp) != nr) {
			fprintf(stderr,"error reading from file\n");
			exit(0);
		}
		count -= nr;
		buf += nr;
	}
	fclose(fp);
	return buf0;
}

static void bufclr(s,n)
char *s;
long n;
{
	while (n-- > 0) *s++ = 0;
}

char *ralloc(size)
long size;
{
char *ptr;
	if ((ptr = malloc(size)) == NULL) {
		fprintf(stderr,"can't malloc %ld bytes\n",size);
		exit(0);
	}
	bufclr(ptr,size);
	return ptr;
}
