/*
 * Author: Bill Gardner
 * Copyright (c) 1995 MIT Media Lab All Rights Reserved.
 */
void s_split(short *stereo, long ns, short *left, short *right);
void c_opt_mul(complex *a, complex *b, complex *c, long n);
short *file_load(char *fname, long *p_n);
char *ralloc(long size);


