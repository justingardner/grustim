/*
 * Basic array manipulation: clear, copy, and allocate.
 *
 * Author: Bill Gardner
 * Copyright (c) 1995 MIT Media Lab All Rights Reserved.
 */
#define _ARRAYLIB_H_

extern void b_clear(char *, long);
extern void s_clear(short *, long);
extern void l_clear(long *, long);
extern void f_clear(float *, long);
extern void c_clear(complex *, long);

extern void b_copy(char *, char *, long);
extern void s_copy(short *, short *, long);
extern void l_copy(long *, long *, long);
extern void f_copy(float *, float *, long);
extern void c_copy(complex *, complex *, long);

extern void l_s_copy(long *, short *, long);
extern void s_l_copy(short *, long *, long);
extern void l_f_copy(long *, float *, long);
extern void f_l_copy(float *, long *, long);
extern void f_s_copy(float *, short *, long);
extern void s_f_copy(short *, float *, long);

extern char *b_alloc(long);
extern short *s_alloc(long);
extern long *l_alloc(long);
extern float *f_alloc(long);
extern complex *c_alloc(long);

extern int s_round(double);
extern long l_round(double);
