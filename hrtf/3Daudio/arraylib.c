/*
 * Basic array manipulation: clear, copy, and allocate.
 *
 * Author: Bill Gardner
 * Copyright (c) 1995 MIT Media Lab All Rights Reserved.
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	"complex.h"
#include	"arraylib.h"

/*
 * Clearing.
 */
static long *_l_clear(p, n)
register long *p;
register long n;
{
	while (n-- > 0) *p++ = 0;
	return p;
}

static void _b_clear(p, n)
char *p;
long n;
{
	while (n-- > 0) *p++ = 0;
}

/*
 * bufclr optimized for big arrays.
 */
void b_clear(p, n)
char *p;
long n;
{
	p = (char *) _l_clear((long *) p, n >> 2);
	_b_clear(p, n & 3);
}

void s_clear(p, n)
short *p;
long n;
{
	b_clear((char *) p, (long) n * sizeof(short));
}

void l_clear(p, n)
long *p;
long n;
{
	b_clear((char *) p, (long)  n * sizeof(long));
}

void f_clear(p, n)
float *p;
long n;
{
	b_clear((char *) p, (long)  n * sizeof(float));
}

void c_clear(p, n)
complex *p;
long n;
{
	b_clear((char *) p, (long)  n * sizeof(complex));
}

/*
 * Copying.
 */

/*
 * This macro for copying deals with the case of overlapped
 * copies to a higher address which must be done backwards
 * from the ends of the arrays.
 */
#define COPY_MACRO(t, f, n)	\
	if ((char *) t > (char *) f) {		\
		t += n;			\
		f += n;			\
		while (n-- > 0) *--t = *--f;	\
	}					\
	else				\
		while (n-- > 0) *t++ = *f++;

void b_copy(t, f, n)
char *t, *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void s_copy(t, f, n)
short *t, *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void l_copy(t, f, n)
long *t, *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void f_copy(t, f, n)
float *t, *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void c_copy(t, f, n)
complex *t, *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void s_l_copy(t, f, n)
short *t;
long *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void l_s_copy(t, f, n)
long *t;
short *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void f_l_copy(t, f, n)
float *t;
long *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void l_f_copy(t, f, n)
long *t;
float *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void f_s_copy(t, f, n)
float *t;
short *f;
long n;
{
	COPY_MACRO(t, f, n);
}

void s_f_copy(t, f, n)
short *t;
float *f;
long n;
{
	if ((char *) t > (char *) f) {
		t += n;
		f += n;
		while (n-- > 0) *--t = s_round(*--f);
	}
	else
		while (n-- > 0) *t++ = s_round(*f++);
}

/*
 * Allocating.
 */
char *b_alloc(n)
long n;
{
char *p;
	if ((p = malloc(n)) == 0L) {
		fprintf(stderr,"can't allocate %ld bytes\n",n);
		exit(0);
	}
	b_clear(p,n);
	return p;
}

short *s_alloc(n)
long n;
{
	return (short *) b_alloc(n * sizeof(short));
}

long *l_alloc(n)
long n;
{
	return (long *) b_alloc(n * sizeof(long));
}

float *f_alloc(n)
long n;
{
	return (float *) b_alloc(n * sizeof(float));
}

complex *c_alloc(n)
long n;
{
	return (complex *) b_alloc(n * sizeof(complex));
}

/*
 * Rounding. Needed by copiers, that's why it's here.
 */
int s_round(f)
double f;
{
	return (f > 0) ? f + 0.5 : f - 0.5;
}

long l_round(f)
double f;
{
	return (f > 0) ? f + 0.5 : f - 0.5;
}
