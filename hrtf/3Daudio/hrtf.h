/*
 * Author: Bill Gardner
 * Copyright (c) 1995 MIT Media Lab All Rights Reserved.
 */

void get_hrtf(double elev, double azim, complex **p_left, complex **p_right);
void free_hrtfs();
void do_ctl(int ctl_num, int ctl_val);
void get_cur_hrtf(complex **p_left, complex **p_right, double *p_gain);
void do_init(float **p_ibuf);
