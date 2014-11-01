/*
 * rtl-sdr, turns your Realtek RTL2832 based DVB dongle into a SDR receiver
 * Copyright (C) 2012 by Steve Markgraf <steve@steve-m.de>
 * Copyright (C) 2012 by Hoernchen <la@tfc-server.de>
 * Copyright (C) 2012 by Kyle Keen <keenerd@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/*
 * rtl_power: general purpose FFT integrator
 * -f low_freq:high_freq:max_bin_size
 * -i seconds
 * outputs CSV
 * time, low, high, step, db, db, db ...
 * db optional?  raw output might be better for noise correction
 * todo:
 *	threading
 *	randomized hopping
 *	noise correction
 *	continuous IIR
 *	general astronomy usefulness
 *	multiple dongles
 *	multiple FFT workers
 *	check edge cropping for off-by-one and rounding errors
 *	1.8MS/s for hiding xtal harmonics
 */

#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>
#include <fcntl.h>
#include <io.h>
#include "getopt/getopt.h"
#define usleep(x) Sleep(x/1000)
#if defined(_MSC_VER) && _MSC_VER < 1800
#define round(x) (x > 0.0 ? floor(x + 0.5): ceil(x - 0.5))
#endif
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <pthread.h>
#include <libusb.h>

#include "rtl-sdr.h"
#include "convenience/convenience.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define DEFAULT_BUF_LENGTH		(1 * 16384)
#define AUTO_GAIN			-100
#define BUFFER_DUMP			(1<<12)

#define MAXIMUM_RATE			3200000
#define MINIMUM_RATE			1000000
#define DEFAULT_TARGET			2400000

#define MAXIMUM_FFT			32

static volatile int do_exit = 0;
static rtlsdr_dev_t *dev = NULL;
FILE *file;

struct sine_table
{
	int16_t* Sinewave;
	int N_WAVE;
	int LOG2_N_WAVE;
};

struct sine_table s_tables[MAXIMUM_FFT];

struct tuning_state
/* one per tuning range */
{
	int freq;
	int rate;
	int gain;
	int bin_e;
	int16_t *fft_buf;
	int64_t *avg;  /* length == 2^bin_e */
	int samples;
	int downsample;
	int downsample_passes;  /* for the recursive filter */
	int comp_fir_size;
	int peak_hold;
	int linear;
	int bin_spec;
	double crop;
	int crop_i1, crop_i2;
	int freq_low, freq_high;
	//pthread_rwlock_t avg_lock;
	//pthread_mutex_t avg_mutex;
	/* having the iq buffer here is wasteful, but will avoid contention */
	uint8_t *buf8;
	int buf_len;
	//int *comp_fir;
	//pthread_rwlock_t buf_lock;
	//pthread_mutex_t buf_mutex;
	int *window_coefs;
	struct sine_table *sine;  /* points to an element of s_tables */
};

struct channel_solve
/* details required to find optimal tuning */
{
	int upper, lower, bin_spec;
	int hops, bw_wanted, bw_needed;
	int bin_e, downsample, downsample_passes;
	double crop, crop_tmp;
};

struct misc_settings
{
	int boxcar;
	int comp_fir_size;
	int peak_hold;
	int linear;
	int target_rate;
	double crop;
	int gain;
	double (*window_fn)(int, int);
	int smoothing;
};

/* 3000 is enough for 3GHz b/w worst case */
#define MAX_TUNES	4000
struct tuning_state tunes[MAX_TUNES];
int tune_count = 0;

void usage(void)
{
	fprintf(stderr,
		"rtl_power, a simple FFT logger for RTL2832 based DVB-T receivers\n\n"
		"Use:\trtl_power -f freq_range [-options] [-f freq2 -opts2] [filename]\n"
		"\t-f lower:upper:bin_size [Hz]\n"
		"\t  valid range for bin_size is 1Hz - 2.8MHz\n"
		"\t  multiple frequency ranges are supported\n"
		"\t[-i integration_interval (default: 10 seconds)]\n"
		"\t  buggy if a full sweep takes longer than the interval\n"
		"\t[-1 enables single-shot mode (default: off)]\n"
		"\t[-e exit_timer (default: off/0)]\n"
		//"\t[-s avg/iir smoothing (default: avg)]\n"
		//"\t[-t threads (default: 1)]\n"
		"\t[-d device_index (default: 0)]\n"
		"\t[-g tuner_gain (default: automatic)]\n"
		"\t[-p ppm_error (default: 0)]\n"
		"\tfilename (a '-' dumps samples to stdout)\n"
		"\t  omitting the filename also uses stdout\n"
		"\n"
		"Experimental options:\n"
		"\t[-w window (default: rectangle)]\n"
		"\t  hamming, blackman, blackman-harris, hann-poisson, bartlett, youssef\n"
		// kaiser
		"\t[-c crop_percent (default: 0%% suggested: 20%%)]\n"
		"\t  discards data at the edges, 100%% discards everything\n"
		"\t  has no effect for bins larger than 1MHz\n"
		"\t  this value is a minimum crop size, more may be discarded\n"
		"\t[-F fir_size (default: disabled)]\n"
		"\t  enables low-leakage downsample filter,\n"
		"\t  fir_size can be 0 or 9.  0 has bad roll off,\n"
		"\t  try -F 0 with '-c 50%%' to hide the roll off\n"
		"\t[-r max_sample_rate (default: 2.4M)]\n"
		"\t  possible values are 2M to 3.2M\n"
		"\t[-P enables peak hold (default: off/averaging)]\n"
		"\t[-L enable linear output (default: off/dB)]\n"
		"\t[-D direct_sampling_mode, 0 (default/off), 1 (I), 2 (Q), 3 (no-mod)]\n"
		"\t[-O enable offset tuning (default: off)]\n"
		"\n"
		"CSV FFT output columns:\n"
		"\tdate, time, Hz low, Hz high, Hz step, samples, dbm, dbm, ...\n\n"
		"Examples:\n"
		"\trtl_power -f 88M:108M:125k fm_stations.csv\n"
		"\t  creates 160 bins across the FM band,\n"
		"\t  individual stations should be visible\n"
		"\trtl_power -f 100M:1G:1M -i 5m -1 survey.csv\n"
		"\t  a five minute low res scan of nearly everything\n"
		"\trtl_power -f ... -i 15m -1 log.csv\n"
		"\t  integrate for 15 minutes and exit afterwards\n"
		"\trtl_power -f ... -e 1h | gzip > log.csv.gz\n"
		"\t  collect data for one hour and compress it on the fly\n\n"
		"\tIf you have issues writing +2GB logs on a 32bit platform\n"
		"\tuse redirection (rtl_power ... > filename.csv) instead\n\n"
		"Convert CSV to a waterfall graphic with:\n"
		"  https://github.com/keenerd/rtl-sdr-misc/blob/master/heatmap/heatmap.py \n"
		"More examples at http://kmkeen.com/rtl-power/\n");
	exit(1);
}

void multi_bail(void)
{
	if (do_exit == 1)
	{
		fprintf(stderr, "Signal caught, finishing scan pass.\n");
	}
	if (do_exit >= 2)
	{
		fprintf(stderr, "Signal caught, aborting immediately.\n");
	}
}

#ifdef _WIN32
BOOL WINAPI
sighandler(int signum)
{
	if (CTRL_C_EVENT == signum) {
		do_exit++;
		multi_bail();
		return TRUE;
	}
	return FALSE;
}
#else
static void sighandler(int signum)
{
	do_exit++;
	multi_bail();
}
#endif

/* more cond dumbness */
#define safe_cond_signal(n, m) pthread_mutex_lock(m); pthread_cond_signal(n); pthread_mutex_unlock(m)
#define safe_cond_wait(n, m) pthread_mutex_lock(m); pthread_cond_wait(n, m); pthread_mutex_unlock(m)

/* {length, coef, coef, coef}  and scaled by 2^15
   for now, only length 9, optimal way to get +85% bandwidth */
#define CIC_TABLE_MAX 10
int cic_9_tables[][10] = {
	{0,},
	{9, -156,  -97, 2798, -15489, 61019, -15489, 2798,  -97, -156},
	{9, -128, -568, 5593, -24125, 74126, -24125, 5593, -568, -128},
	{9, -129, -639, 6187, -26281, 77511, -26281, 6187, -639, -129},
	{9, -122, -612, 6082, -26353, 77818, -26353, 6082, -612, -122},
	{9, -120, -602, 6015, -26269, 77757, -26269, 6015, -602, -120},
	{9, -120, -582, 5951, -26128, 77542, -26128, 5951, -582, -120},
	{9, -119, -580, 5931, -26094, 77505, -26094, 5931, -580, -119},
	{9, -119, -578, 5921, -26077, 77484, -26077, 5921, -578, -119},
	{9, -119, -577, 5917, -26067, 77473, -26067, 5917, -577, -119},
	{9, -199, -362, 5303, -25505, 77489, -25505, 5303, -362, -199},
};

#if defined(_MSC_VER) && _MSC_VER < 1800
double log2(double n)
{
	return log(n) / log(2.0);
}
#endif

/* FFT based on fix_fft.c by Roberts, Slaney and Bouras
   http://www.jjj.de/fft/fftpage.html
   16 bit ints for everything
   -32768..+32768 maps to -1.0..+1.0
*/

void sine_table(int size)
{
	int i;
	double d;
	struct sine_table *sine;
	if (size > (MAXIMUM_FFT-1)) {
		fprintf(stderr, "Maximum FFT is 2^%i\n", MAXIMUM_FFT-1);
		exit(1);
	}
	sine = &s_tables[size];
	if (sine->LOG2_N_WAVE == size) {
		return;}
	sine->LOG2_N_WAVE = size;
	sine->N_WAVE = 1 << sine->LOG2_N_WAVE;
	sine->Sinewave = malloc(sizeof(int16_t) * sine->N_WAVE*3/4);
	for (i=0; i<sine->N_WAVE*3/4; i++)
	{
		d = (double)i * 2.0 * M_PI / sine->N_WAVE;
		sine->Sinewave[i] = (int)round(32767*sin(d));
		//printf("%i\n", sine->Sinewave[i]);
	}
}

void generate_sine_tables(void)
{
	struct tuning_state *ts;
	int i;
	for (i=0; i < tune_count; i++) {
		ts = &tunes[i];
		sine_table(ts->bin_e);
		ts->sine = &s_tables[ts->bin_e];
		ts->fft_buf = malloc(ts->buf_len * sizeof(int16_t));
	}
}

inline int16_t FIX_MPY(int16_t a, int16_t b)
/* fixed point multiply and scale */
{
	int c = ((int)a * (int)b) >> 14;
	b = c & 0x01;
	return (c >> 1) + b;
}

int fix_fft(int16_t iq[], int m, struct sine_table *sine)
/* interleaved iq[], 0 <= n < 2**m, changes in place */
{
	int mr, nn, i, j, l, k, istep, n, shift;
	int16_t qr, qi, tr, ti, wr, wi;
	n = 1 << m;
	if (n > sine->N_WAVE)
		{return -1;}
	mr = 0;
	nn = n - 1;
	/* decimation in time - re-order data */
	for (m=1; m<=nn; ++m) {
		l = n;
		do
			{l >>= 1;}
		while (mr+l > nn);
		mr = (mr & (l-1)) + l;
		if (mr <= m)
			{continue;}
		// real = 2*m, imag = 2*m+1
		tr = iq[2*m];
		iq[2*m] = iq[2*mr];
		iq[2*mr] = tr;
		ti = iq[2*m+1];
		iq[2*m+1] = iq[2*mr+1];
		iq[2*mr+1] = ti;
	}
	l = 1;
	k = sine->LOG2_N_WAVE-1;
	while (l < n) {
		shift = 1;
		istep = l << 1;
		for (m=0; m<l; ++m) {
			j = m << k;
			wr =  sine->Sinewave[j+sine->N_WAVE/4];
			wi = -sine->Sinewave[j];
			if (shift) {
				wr >>= 1; wi >>= 1;}
			for (i=m; i<n; i+=istep) {
				j = i + l;
				tr = FIX_MPY(wr,iq[2*j]) - FIX_MPY(wi,iq[2*j+1]);
				ti = FIX_MPY(wr,iq[2*j+1]) + FIX_MPY(wi,iq[2*j]);
				qr = iq[2*i];
				qi = iq[2*i+1];
				if (shift) {
					qr >>= 1; qi >>= 1;}
				iq[2*j] = qr - tr;
				iq[2*j+1] = qi - ti;
				iq[2*i] = qr + tr;
				iq[2*i+1] = qi + ti;
			}
		}
		--k;
		l = istep;
	}
	return 0;
}

double rectangle(int i, int length)
{
	return 1.0;
}

double hamming(int i, int length)
{
	double a, b, w, N1;
	a = 25.0/46.0;
	b = 21.0/46.0;
	N1 = (double)(length-1);
	w = a - b*cos(2*i*M_PI/N1);
	return w;
}

double blackman(int i, int length)
{
	double a0, a1, a2, w, N1;
	a0 = 7938.0/18608.0;
	a1 = 9240.0/18608.0;
	a2 = 1430.0/18608.0;
	N1 = (double)(length-1);
	w = a0 - a1*cos(2*i*M_PI/N1) + a2*cos(4*i*M_PI/N1);
	return w;
}

double blackman_harris(int i, int length)
{
	double a0, a1, a2, a3, w, N1;
	a0 = 0.35875;
	a1 = 0.48829;
	a2 = 0.14128;
	a3 = 0.01168;
	N1 = (double)(length-1);
	w = a0 - a1*cos(2*i*M_PI/N1) + a2*cos(4*i*M_PI/N1) - a3*cos(6*i*M_PI/N1);
	return w;
}

double hann_poisson(int i, int length)
{
	double a, N1, w;
	a = 2.0;
	N1 = (double)(length-1);
	w = 0.5 * (1 - cos(2*M_PI*i/N1)) * \
	    pow(M_E, (-a*(double)abs((int)(N1-1-2*i)))/N1);
	return w;
}

double youssef(int i, int length)
/* really a blackman-harris-poisson window, but that is a mouthful */
{
	double a, a0, a1, a2, a3, w, N1;
	a0 = 0.35875;
	a1 = 0.48829;
	a2 = 0.14128;
	a3 = 0.01168;
	N1 = (double)(length-1);
	w = a0 - a1*cos(2*i*M_PI/N1) + a2*cos(4*i*M_PI/N1) - a3*cos(6*i*M_PI/N1);
	a = 0.0025;
	w *= pow(M_E, (-a*(double)abs((int)(N1-1-2*i)))/N1);
	return w;
}

double kaiser(int i, int length)
// todo, become more smart
{
	return 1.0;
}

double bartlett(int i, int length)
{
	double N1, L, w;
	L = (double)length;
	N1 = L - 1;
	w = (i - N1/2) / (L/2);
	if (w < 0) {
		w = -w;}
	w = 1 - w;
	return w;
}

void rms_power(struct tuning_state *ts)
/* for bins between 1MHz and 2MHz */
{
	int i, s;
	uint8_t *buf = ts->buf8;
	int buf_len = ts->buf_len;
	int64_t p, t;
	double dc, err;

	p = t = 0L;
	for (i=0; i<buf_len; i++) {
		s = (int)buf[i] - 127;
		t += (int64_t)s;
		p += (int64_t)(s * s);
	}
	/* correct for dc offset in squares */
	dc = (double)t / (double)buf_len;
	err = t * 2 * dc - dc * dc * buf_len;
	p -= (int64_t)round(err);

	if (!ts->peak_hold) {
		ts->avg[0] += p;
	} else {
		ts->avg[0] = MAX(ts->avg[0], p);
	}
	ts->samples += 1;
}

/* todo, add errors to parse_freq, solve_foo */

int parse_frequency(char *arg, struct channel_solve *c)
{
	char *start, *stop, *step;
	/* hacky string parsing */
	start = arg;
	stop = strchr(start, ':') + 1;
	stop[-1] = '\0';
	step = strchr(stop, ':') + 1;
	step[-1] = '\0';
	c->lower = (int)atofs(start);
	c->upper = (int)atofs(stop);
	c->bin_spec = (int)atofs(step);
	stop[-1] = ':';
	step[-1] = ':';
	return 0;
}

int solve_giant_bins(struct channel_solve *c)
{
	c->bw_wanted = c->bin_spec;
	c->bw_needed = c->bin_spec;
	c->hops = (c->upper - c->lower) / c->bin_spec;
	c->bin_e = 0;
	c->crop_tmp = 0;
	return 0;
}

int solve_downsample(struct channel_solve *c, int target_rate, int boxcar)
{
	int scan_size, bins_wanted, bins_needed, ds_next, bw;

	scan_size = c->upper - c->lower;
	c->hops = 1;
	c->bw_wanted = scan_size;

	bins_wanted = (int)ceil((double)scan_size / (double)c->bin_spec);
	c->bin_e = (int)ceil(log2(bins_wanted));
	while (1) {
		bins_needed = 1 << c->bin_e;
		c->crop_tmp = (double)(bins_needed - bins_wanted) / (double)bins_needed;
		if (c->crop_tmp >= c->crop) {
			break;}
		c->bin_e++;
	}

	c->downsample = 1;
	c->downsample_passes = 0;
	while (1) {
		bw = (int)((double)scan_size / (1.0 - c->crop_tmp));
		c->bw_needed = bw * c->downsample;

		if (boxcar) {
			ds_next = c->downsample + 1;
		} else {
			ds_next = c->downsample * 2;
		}
		if ((bw * ds_next) > target_rate) {
			break;}

		c->downsample = ds_next;
		if (!boxcar) {
			c->downsample_passes++;}
	}

	return 0;
}

int solve_hopping(struct channel_solve *c, int target_rate)
{
	int i, scan_size, bins_all, bins_crop, bins_2, min_hops;
	scan_size = c->upper - c->lower;
	min_hops = scan_size / MAXIMUM_RATE - 1;
	if (min_hops < 1) {
		min_hops = 1;}
	/* evenly sized ranges, as close to target_rate as possible */
	for (i=min_hops; i<MAX_TUNES; i++) {
		c->bw_wanted = scan_size / i;
		bins_all = scan_size / c->bin_spec;
		bins_crop = (int)ceil((double)bins_all / (double)i);
		c->bin_e = (int)ceil(log2(bins_crop));
		bins_2 = 1 << c->bin_e;
		c->bw_needed = bins_2 * c->bin_spec;
		c->crop_tmp = (double)(bins_2 - bins_crop) / (double)bins_2;
		if (c->bw_needed > target_rate) {
			continue;}
		if (c->crop_tmp < c->crop) {
			continue;}
		c->hops = i;
		break;
	}
	return 0;
}

void frequency_range(char *arg, struct misc_settings *ms)
/* flesh out the tunes[] for scanning */
{
	struct channel_solve c;
	struct tuning_state *ts;
	int i, j, buf_len, length, hop_bins, logged_bins, planned_bins;
	int lower_edge, actual_bw, upper_perfect, remainder;

	fprintf(stderr, "Range: %s\n", arg);
	parse_frequency(arg, &c);
	c.downsample = 1;
	c.downsample_passes = 0;
	c.crop = ms->crop;

	if (ms->target_rate < 2 * MINIMUM_RATE) {
		ms->target_rate = 2 * MINIMUM_RATE;
	}
	if (ms->target_rate > MAXIMUM_RATE) {
		ms->target_rate = MAXIMUM_RATE;
	}
	if ((ms->crop < 0.0) || (ms->crop > 1.0)) {
		fprintf(stderr, "Crop value outside of 0 to 1.\n");
		exit(1);
	}


	if (c.bin_spec >= MINIMUM_RATE) {
		fprintf(stderr, "Mode: rms power\n");
		solve_giant_bins(&c);
	} else if ((c.upper - c.lower) < MINIMUM_RATE) {
		fprintf(stderr, "Mode: downsampling\n");
		solve_downsample(&c, ms->target_rate, ms->boxcar);
	} else {
		fprintf(stderr, "Mode: normal\n");
		solve_hopping(&c, ms->target_rate);
	}
	c.crop = c.crop_tmp;

	if ((tune_count+c.hops) > MAX_TUNES) {
		fprintf(stderr, "Error: bandwidth too wide.\n");
		exit(1);
	}
	buf_len = 2 * (1<<c.bin_e) * c.downsample;
	if (buf_len < DEFAULT_BUF_LENGTH) {
		buf_len = DEFAULT_BUF_LENGTH;
	}
	/* build the array */
	logged_bins = 0;
	lower_edge = c.lower;
	planned_bins = (c.upper - c.lower) / c.bin_spec;
	for (i=0; i < c.hops; i++) {
		ts = &tunes[tune_count + i];
		/* copy common values */
		ts->rate = c.bw_needed;
		ts->gain = ms->gain;
		ts->bin_e = c.bin_e;
		ts->samples = 0;
		ts->bin_spec = c.bin_spec;
		ts->crop = c.crop;
		ts->downsample = c.downsample;
		ts->downsample_passes = c.downsample_passes;
		ts->comp_fir_size = ms->comp_fir_size;
		ts->peak_hold = ms->peak_hold;
		ts->linear = ms->linear;
		ts->avg = (int64_t*)malloc((1<<c.bin_e) * sizeof(int64_t));
		if (!ts->avg) {
			fprintf(stderr, "Error: malloc.\n");
			exit(1);
		}
		for (j=0; j<(1<<c.bin_e); j++) {
			ts->avg[j] = 0L;
		}
		ts->buf8 = (uint8_t*)malloc(buf_len * sizeof(uint8_t));
		if (!ts->buf8) {
			fprintf(stderr, "Error: malloc.\n");
			exit(1);
		}
		ts->buf_len = buf_len;
		length = 1 << c.bin_e;
		ts->window_coefs = malloc(length * sizeof(int));
		for (j=0; j<length; j++) {
			ts->window_coefs[j] = (int)(256*ms->window_fn(j, length));
		}
		/* calculate unique values */
		ts->freq_low = lower_edge;
		hop_bins = c.bw_wanted / c.bin_spec;
		actual_bw = hop_bins * c.bin_spec;
		ts->freq_high = lower_edge + actual_bw;
		upper_perfect = c.lower + (i+1) * c.bw_wanted;
		if (ts->freq_high + c.bin_spec <= upper_perfect) {
			hop_bins += 1;
			actual_bw = hop_bins * c.bin_spec;
			ts->freq_high = lower_edge + actual_bw;
		}
		remainder = planned_bins - logged_bins - hop_bins;
		if (i == c.hops-1 && remainder > 0) {
			hop_bins += remainder;
			actual_bw = hop_bins * c.bin_spec;
			ts->freq_high = lower_edge + actual_bw;
		}
		logged_bins += hop_bins;
		ts->crop_i1 = (length - hop_bins) / 2;
		ts->crop_i2 = ts->crop_i1 + hop_bins - 1;
		ts->freq = (lower_edge - ts->crop_i1 * c.bin_spec) + c.bw_needed/(2*c.downsample);
		/* prep for next hop */
		lower_edge = ts->freq_high;
	}
	tune_count += c.hops;
	/* report */
	fprintf(stderr, "Number of frequency hops: %i\n", c.hops);
	fprintf(stderr, "Dongle bandwidth: %iHz\n", c.bw_needed);
	fprintf(stderr, "Downsampling by: %ix\n", c.downsample);
	fprintf(stderr, "Cropping by: %0.2f%%\n", c.crop*100);
	fprintf(stderr, "Total FFT bins: %i\n", c.hops * (1<<c.bin_e));
	fprintf(stderr, "Logged FFT bins: %i\n", logged_bins);
	fprintf(stderr, "FFT bin size: %iHz\n", c.bin_spec);
	fprintf(stderr, "Buffer size: %i bytes (%0.2fms)\n", buf_len, 1000 * 0.5 * (float)buf_len / (float)c.bw_needed);
	fprintf(stderr, "\n");
}

void retune(rtlsdr_dev_t *d, int freq)
{
	uint8_t dump[BUFFER_DUMP];
	int f, n_read;
	f = (int)rtlsdr_get_center_freq(d);
	if (f == freq) {
		return;}
	rtlsdr_set_center_freq(d, (uint32_t)freq);
	/* wait for settling and flush buffer */
	usleep(5000);
	rtlsdr_read_sync(d, &dump, BUFFER_DUMP, &n_read);
	if (n_read != BUFFER_DUMP) {
		fprintf(stderr, "Error: bad retune.\n");}
}

void rerate(rtlsdr_dev_t *d, int rate)
{
	uint32_t r;
	r = rtlsdr_get_sample_rate(d);
	if ((int)r == rate) {
		return;}
	rtlsdr_set_sample_rate(dev, (uint32_t)rate);
}

void regain(rtlsdr_dev_t *d, int gain)
{
	int g;
	g = rtlsdr_get_tuner_gain(dev);
	if (g == gain) {
		return;}
	if (gain == AUTO_GAIN) {
		/* switch to auto */
		rtlsdr_set_tuner_gain_mode(dev, 0);
		return;
	}
	if (g == AUTO_GAIN) {
		/* switch to manual */
		rtlsdr_set_tuner_gain_mode(dev, 1);
	}
	rtlsdr_set_tuner_gain(dev, gain);
}

void fifth_order(int16_t *data, int length)
/* for half of interleaved data */
{
	int i;
	int a, b, c, d, e, f;
	a = data[0];
	b = data[2];
	c = data[4];
	d = data[6];
	e = data[8];
	f = data[10];
	/* a downsample should improve resolution, so don't fully shift */
	/* ease in instead of being stateful */
	data[0] = ((a+b)*10 + (c+d)*5 + d + f) >> 4;
	data[2] = ((b+c)*10 + (a+d)*5 + e + f) >> 4;
	data[4] = (a + (b+e)*5 + (c+d)*10 + f) >> 4;
	for (i=12; i<length; i+=4) {
		a = c;
		b = d;
		c = e;
		d = f;
		e = data[i-2];
		f = data[i];
		data[i/2] = (a + (b+e)*5 + (c+d)*10 + f) >> 4;
	}
}

void remove_dc(int16_t *data, int length)
/* works on interleaved data */
{
	int i;
	int16_t ave;
	int64_t sum = 0L;
	for (i=0; i < length; i+=2) {
		sum += data[i];
	}
	ave = (int16_t)(sum / (int64_t)(length));
	if (ave == 0) {
		return;}
	for (i=0; i < length; i+=2) {
		data[i] -= ave;
	}
}

void generic_fir(int16_t *data, int length, int *fir)
/* Okay, not at all generic.  Assumes length 9, fix that eventually. */
{
	int d, temp, sum;
	int hist[9] = {0,};
	/* cheat on the beginning, let it go unfiltered */
	for (d=0; d<18; d+=2) {
		hist[d/2] = data[d];
	}
	for (d=18; d<length; d+=2) {
		temp = data[d];
		sum = 0;
		sum += (hist[0] + hist[8]) * fir[1];
		sum += (hist[1] + hist[7]) * fir[2];
		sum += (hist[2] + hist[6]) * fir[3];
		sum += (hist[3] + hist[5]) * fir[4];
		sum +=            hist[4]  * fir[5];
		data[d] = (int16_t)(sum >> 15) ;
		hist[0] = hist[1];
		hist[1] = hist[2];
		hist[2] = hist[3];
		hist[3] = hist[4];
		hist[4] = hist[5];
		hist[5] = hist[6];
		hist[6] = hist[7];
		hist[7] = hist[8];
		hist[8] = temp;
	}
}

void downsample_iq(int16_t *data, int length)
{
	fifth_order(data, length);
	//remove_dc(data, length);
	fifth_order(data+1, length-1);
	//remove_dc(data+1, length-1);
}

int64_t real_conj(int16_t real, int16_t imag)
/* real(n * conj(n)) */
{
	return ((int64_t)real*(int64_t)real + (int64_t)imag*(int64_t)imag);
}

void scanner(void)
{
	int i, j, j2, n_read, offset, bin_e, bin_len, buf_len, ds, ds_p;
	int32_t w;
	struct tuning_state *ts;
	int16_t *fft_buf;
	bin_e = tunes[0].bin_e;
	bin_len = 1 << bin_e;
	buf_len = tunes[0].buf_len;
	for (i=0; i<tune_count; i++) {
		if (do_exit >= 2)
			{return;}
		ts = &tunes[i];
		fft_buf = ts->fft_buf;
		regain(dev, ts->gain);
		rerate(dev, ts->rate);
		retune(dev, ts->freq);
		rtlsdr_read_sync(dev, ts->buf8, buf_len, &n_read);
		if (n_read != buf_len) {
			fprintf(stderr, "Error: dropped samples.\n");}
		/* rms */
		if (bin_len == 1) {
			rms_power(ts);
			continue;
		}
		/* prep for fft */
		for (j=0; j<buf_len; j++) {
			fft_buf[j] = (int16_t)ts->buf8[j] - 127;
		}
		ds = ts->downsample;
		ds_p = ts->downsample_passes;
		if (ds_p) {  /* recursive */
			for (j=0; j < ds_p; j++) {
				downsample_iq(fft_buf, buf_len >> j);
			}
			/* droop compensation */
			if (ts->comp_fir_size == 9 && ds_p <= CIC_TABLE_MAX) {
				generic_fir(fft_buf, buf_len >> j, cic_9_tables[ds_p]);
				generic_fir(fft_buf+1, (buf_len >> j)-1, cic_9_tables[ds_p]);
			}
		} else if (ds > 1) {  /* boxcar */
			j=2, j2=0;
			while (j < buf_len) {
				fft_buf[j2]   += fft_buf[j];
				fft_buf[j2+1] += fft_buf[j+1];
				fft_buf[j] = 0;
				fft_buf[j+1] = 0;
				j += 2;
				if (j % (ds*2) == 0) {
					j2 += 2;}
			}
		}
		remove_dc(fft_buf, buf_len / ds);
		remove_dc(fft_buf+1, (buf_len / ds) - 1);
		/* window function and fft */
		for (offset=0; offset<(buf_len/ds); offset+=(2*bin_len)) {
			// todo, let rect skip this
			for (j=0; j<bin_len; j++) {
				w =  (int32_t)fft_buf[offset+j*2];
				w *= (int32_t)(ts->window_coefs[j]);
				//w /= (int32_t)(ds);
				fft_buf[offset+j*2]   = (int16_t)w;
				w =  (int32_t)fft_buf[offset+j*2+1];
				w *= (int32_t)(ts->window_coefs[j]);
				//w /= (int32_t)(ds);
				fft_buf[offset+j*2+1] = (int16_t)w;
			}
			fix_fft(fft_buf+offset, bin_e, ts->sine);
			if (!ts->peak_hold) {
				for (j=0; j<bin_len; j++) {
					ts->avg[j] += real_conj(fft_buf[offset+j*2], fft_buf[offset+j*2+1]);
				}
			} else {
				for (j=0; j<bin_len; j++) {
					ts->avg[j] = MAX(real_conj(fft_buf[offset+j*2], fft_buf[offset+j*2+1]), ts->avg[j]);
				}
			}
			ts->samples += ds;
		}
	}
}

void csv_dbm(struct tuning_state *ts)
{
	int i, len;
	int64_t tmp;
	double dbm;
	char *sep = ", ";
	len = 1 << ts->bin_e;
	/* fix FFT stuff quirks */
	if (ts->bin_e > 0) {
		/* nuke DC component (not effective for all windows) */
		ts->avg[0] = ts->avg[1];
		/* FFT is translated by 180 degrees */
		for (i=0; i<len/2; i++) {
			tmp = ts->avg[i];
			ts->avg[i] = ts->avg[i+len/2];
			ts->avg[i+len/2] = tmp;
		}
	}
	/* Hz low, Hz high, Hz step, samples, dbm, dbm, ... */
	fprintf(file, "%i, %i, %i, %i, ", ts->freq_low, ts->freq_high,
		ts->bin_spec, ts->samples);
	// something seems off with the dbm math
	for (i=ts->crop_i1; i<=ts->crop_i2; i++) {
		if (i == ts->crop_i2) {
			sep = "\n";
		}
		dbm  = (double)ts->avg[i];
		dbm /= (double)ts->rate;
		dbm /= (double)ts->samples;
		if (ts->linear) {
			fprintf(file, "%.5g%s", dbm, sep);
		} else {
			dbm = 10 * log10(dbm);
			fprintf(file, "%.2f%s", dbm, sep);
		}
	}
	for (i=0; i<len; i++) {
		ts->avg[i] = 0L;
	}
	ts->samples = 0;
}

void init_misc(struct misc_settings *ms)
{
	ms->target_rate = DEFAULT_TARGET;
	ms->boxcar = 1;
	ms->comp_fir_size = 0;
	ms->crop = 0.0;
	ms->gain = AUTO_GAIN;
	ms->window_fn = rectangle;
	ms->smoothing = 0;
	ms->peak_hold = 0;
	ms->linear = 0;
}

int main(int argc, char **argv)
{
#ifndef _WIN32
	struct sigaction sigact;
#endif
	char *filename = NULL;
	int i, r, opt;
	int f_set = 0;
	int dev_index = 0;
	int dev_given = 0;
	int ppm_error = 0;
	int custom_ppm = 0;
	int interval = 10;
	int fft_threads = 1;
	int single = 0;
	int direct_sampling = 0;
	int offset_tuning = 0;
	char *freq_optarg;
	time_t next_tick;
	time_t time_now;
	time_t exit_time = 0;
	char t_str[50];
	struct tm *cal_time;
	struct misc_settings ms;
	freq_optarg = "";
	init_misc(&ms);

	while ((opt = getopt(argc, argv, "f:i:s:r:t:d:g:p:e:w:c:F:1PLD:Oh")) != -1) {
		switch (opt) {
		case 'f': // lower:upper:bin_size
			if (f_set) {
				frequency_range(freq_optarg, &ms);}
			freq_optarg = strdup(optarg);
			f_set = 1;
			break;
		case 'd':
			dev_index = verbose_device_search(optarg);
			dev_given = 1;
			break;
		case 'g':
			ms.gain = (int)(atof(optarg) * 10);
			break;
		case 'c':
			ms.crop = atofp(optarg);
			break;
		case 'i':
			interval = (int)round(atoft(optarg));
			break;
		case 'e':
			exit_time = (time_t)((int)round(atoft(optarg)));
			break;
		case 's':
			if (strcmp("avg",  optarg) == 0) {
				ms.smoothing = 0;}
			if (strcmp("iir",  optarg) == 0) {
				ms.smoothing = 1;}
			break;
		case 'w':
			if (strcmp("rectangle",  optarg) == 0) {
				ms.window_fn = rectangle;}
			if (strcmp("hamming",  optarg) == 0) {
				ms.window_fn = hamming;}
			if (strcmp("blackman",  optarg) == 0) {
				ms.window_fn = blackman;}
			if (strcmp("blackman-harris",  optarg) == 0) {
				ms.window_fn = blackman_harris;}
			if (strcmp("hann-poisson",  optarg) == 0) {
				ms.window_fn = hann_poisson;}
			if (strcmp("youssef",  optarg) == 0) {
				ms.window_fn = youssef;}
			if (strcmp("kaiser",  optarg) == 0) {
				ms.window_fn = kaiser;}
			if (strcmp("bartlett",  optarg) == 0) {
				ms.window_fn = bartlett;}
			break;
		case 't':
			fft_threads = atoi(optarg);
			break;
		case 'p':
			ppm_error = atoi(optarg);
			custom_ppm = 1;
			break;
		case 'r':
			ms.target_rate = (int)atofs(optarg);
			break;
		case '1':
			single = 1;
			break;
		case 'P':
			ms.peak_hold = 1;
			break;
		case 'L':
			ms.linear = 1;
			break;
		case 'D':
			direct_sampling = atoi(optarg);
			break;
		case 'O':
			offset_tuning = 1;
			break;
		case 'F':
			ms.boxcar = 0;
			ms.comp_fir_size = atoi(optarg);
			break;
		case 'h':
		default:
			usage();
			break;
		}
	}

	if (!f_set) {
		fprintf(stderr, "No frequency range provided.\n");
		exit(1);
	}

	frequency_range(freq_optarg, &ms);

	if (tune_count == 0) {
		usage();}

	if (argc <= optind) {
		filename = "-";
	} else {
		filename = argv[optind];
	}

	if (interval < 1) {
		interval = 1;}

	fprintf(stderr, "Reporting every %i seconds\n", interval);

	if (!dev_given) {
		dev_index = verbose_device_search("0");
	}

	if (dev_index < 0) {
		exit(1);
	}

	r = rtlsdr_open(&dev, (uint32_t)dev_index);
	if (r < 0) {
		fprintf(stderr, "Failed to open rtlsdr device #%d.\n", dev_index);
		exit(1);
	}
#ifndef _WIN32
	sigact.sa_handler = sighandler;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;
	sigaction(SIGINT, &sigact, NULL);
	sigaction(SIGTERM, &sigact, NULL);
	sigaction(SIGQUIT, &sigact, NULL);
	sigaction(SIGPIPE, &sigact, NULL);
	signal(SIGPIPE, SIG_IGN);
#else
	SetConsoleCtrlHandler( (PHANDLER_ROUTINE) sighandler, TRUE );
#endif

	if (direct_sampling) {
		verbose_direct_sampling(dev, direct_sampling);
	}

	if (offset_tuning) {
		verbose_offset_tuning(dev);
	}

	/* Set the tuner gain */
	for (i=0; i<tune_count; i++) {
		if (tunes[i].gain == AUTO_GAIN) {
			continue;}
		tunes[i].gain = nearest_gain(dev, tunes[i].gain);
	}
	if (ms.gain == AUTO_GAIN) {
		verbose_auto_gain(dev);
	} else {
		ms.gain = nearest_gain(dev, ms.gain);
		verbose_gain_set(dev, ms.gain);
	}

	if (!custom_ppm) {
		verbose_ppm_eeprom(dev, &ppm_error);
	}
	verbose_ppm_set(dev, ppm_error);

	if (strcmp(filename, "-") == 0) { /* Write log to stdout */
		file = stdout;
#ifdef _WIN32
		// Is this necessary?  Output is ascii.
		_setmode(_fileno(file), _O_BINARY);
#endif
	} else {
		file = fopen(filename, "wb");
		if (!file) {
			fprintf(stderr, "Failed to open %s\n", filename);
			exit(1);
		}
	}

	generate_sine_tables();

	/* Reset endpoint before we start reading from it (mandatory) */
	verbose_reset_buffer(dev);

	/* actually do stuff */
	rtlsdr_set_sample_rate(dev, (uint32_t)tunes[0].rate);
	next_tick = time(NULL) + interval;
	if (exit_time) {
		exit_time = time(NULL) + exit_time;}
	while (!do_exit) {
		scanner();
		time_now = time(NULL);
		if (time_now < next_tick) {
			continue;}
		// time, Hz low, Hz high, Hz step, samples, dbm, dbm, ...
		cal_time = localtime(&time_now);
		strftime(t_str, 50, "%Y-%m-%d, %H:%M:%S", cal_time);
		for (i=0; i<tune_count; i++) {
			fprintf(file, "%s, ", t_str);
			csv_dbm(&tunes[i]);
		}
		fflush(file);
		while (time(NULL) >= next_tick) {
			next_tick += interval;}
		if (single) {
			do_exit = 1;}
		if (exit_time && time(NULL) >= exit_time) {
			do_exit = 1;}
	}

	/* clean up */

	if (do_exit) {
		fprintf(stderr, "\nUser cancel, exiting...\n");}
	else {
		fprintf(stderr, "\nLibrary error %d, exiting...\n", r);}

	if (file != stdout) {
		fclose(file);}

	rtlsdr_close(dev);
	//free(fft_buf);
	//for (i=0; i<tune_count; i++) {
	//	free(tunes[i].avg);
	//	free(tunes[i].buf8);
	//}
	return r >= 0 ? r : -r;
}

// vim: tabstop=8:softtabstop=8:shiftwidth=8:noexpandtab
