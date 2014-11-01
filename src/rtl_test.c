/*
 * rtl-sdr, turns your Realtek RTL2832 based DVB dongle into a SDR receiver
 * rtl_test, test and benchmark tool
 *
 * Copyright (C) 2012-2014 by Steve Markgraf <steve@steve-m.de>
 * Copyright (C) 2012-2014 by Kyle Keen <keenerd@gmail.com>
 * Copyright (C) 2014 by Michael Tatarinov <kukabu@gmail.com>
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

#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __APPLE__
#include <sys/time.h>
#else
#include <time.h>
#endif

#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>
#include "getopt/getopt.h"
#endif

#include "rtl-sdr.h"
#include "convenience/convenience.h"

#define DEFAULT_SAMPLE_RATE		2048000
#define DEFAULT_BUF_LENGTH		(16 * 16384)
#define MINIMAL_BUF_LENGTH		512
#define MAXIMAL_BUF_LENGTH		(256 * 16384)

#define MHZ(x)	((x)*1000*1000)

#define PPM_DURATION			10
#define PPM_DUMP_TIME			5

#define SCAN_LIMIT			2500000000

static enum {
	NO_BENCHMARK,
	TUNER_BENCHMARK,
	PPM_BENCHMARK
} test_mode = NO_BENCHMARK;

static int do_exit = 0;
static rtlsdr_dev_t *dev = NULL;

static uint32_t samp_rate = DEFAULT_SAMPLE_RATE;

static uint32_t total_samples = 0;
static uint32_t dropped_samples = 0;

static unsigned int ppm_duration = PPM_DURATION;

void usage(void)
{
	fprintf(stderr,
		"rtl_test, a benchmark tool for RTL2832 based DVB-T receivers\n\n"
		"Usage:\n"
		"\t[-s samplerate (default: 2048000 Hz)]\n"
		"\t[-d device_index (default: 0)]\n"
		"\t[-t enable tuner range benchmark]\n"
		"\t[-p[seconds] enable PPM error measurement (default: 10 seconds)]\n"
		"\t[-b output_block_size (default: 16 * 16384)]\n"
		"\t[-S force sync output (default: async)]\n");
	exit(1);
}

#ifdef _WIN32
BOOL WINAPI
sighandler(int signum)
{
	if (CTRL_C_EVENT == signum) {
		fprintf(stderr, "Signal caught, exiting!\n");
		do_exit = 1;
		rtlsdr_cancel_async(dev);
		return TRUE;
	}
	return FALSE;
}
#else
static void sighandler(int signum)
{
	fprintf(stderr, "Signal caught, exiting!\n");
	do_exit = 1;
	rtlsdr_cancel_async(dev);
}
#endif

static void underrun_test(unsigned char *buf, uint32_t len, int mute)
{
	uint32_t i, lost = 0;
	static uint8_t bcnt, uninit = 1;

	if (uninit) {
		bcnt = buf[0];
		uninit = 0;
	}
	for (i = 0; i < len; i++) {
		if(bcnt != buf[i]) {
			lost += (buf[i] > bcnt) ? (buf[i] - bcnt) : (bcnt - buf[i]);
			bcnt = buf[i];
		}

		bcnt++;
	}

	total_samples += len;
	dropped_samples += lost;
	if (mute)
		return;
	if (lost)
		printf("lost at least %d bytes\n", lost);

}

#ifndef _WIN32
static int ppm_gettime(struct timespec *ts)
{
	int rv = ENOSYS;

#ifdef __unix__
	rv = clock_gettime(CLOCK_MONOTONIC, ts);
#elif __APPLE__
	struct timeval tv;

	rv = gettimeofday(&tv, NULL);
	ts->tv_sec = tv.tv_sec;
	ts->tv_nsec = tv.tv_usec * 1000;
#endif
	return rv;
}
#endif

static double get_rate(uint64_t nsamples, int64_t interval, int64_t divider)
{
	return nsamples / (interval / (double)divider);

}

static int ppm_report(uint64_t nsamples, int64_t interval, int64_t divider)
{
	double real_rate = get_rate(nsamples, interval, divider);
	double ppm = 1e6 * (real_rate / (double)samp_rate - 1.);
	return (int)round(ppm);
}

static void ppm_test(uint32_t len)
{
	static uint64_t nsamples = 0;
	static uint64_t nsamples_total = 0;
	static int64_t interval = 0;
	static int64_t interval_total = 0;

#ifdef _WIN32
	static LARGE_INTEGER perfFrequency;
	static LARGE_INTEGER ppm_now;
	static LARGE_INTEGER ppm_recent;
#else
	struct timespec ppm_now;
	static struct timespec ppm_recent;
#endif

	static enum {
		PPM_INIT_NO,
		PPM_INIT_DUMP,
		PPM_INIT_RUN
	} ppm_init = PPM_INIT_NO;

#ifdef _WIN32
	QueryPerformanceCounter(&ppm_now);
#else
	ppm_gettime(&ppm_now);
#endif
	if (ppm_init != PPM_INIT_RUN) {
		/*
		 * Kyle Keen wrote:
		 * PPM_DUMP_TIME throws out the first N seconds of data.
		 * The dongle's PPM is usually very bad when first starting up,
		 * typically incorrect by more than twice the final value.
		 * Discarding the first few seconds allows the value to stabilize much faster.
		*/
		if (ppm_init == PPM_INIT_NO) {
#ifdef _WIN32
			QueryPerformanceFrequency(&perfFrequency);
			ppm_recent.QuadPart = ppm_now.QuadPart + (PPM_DUMP_TIME * perfFrequency.QuadPart);
#else
			ppm_recent.tv_sec = ppm_now.tv_sec + PPM_DUMP_TIME;
#endif
			ppm_init = PPM_INIT_DUMP;
			return;
		}
#ifdef _WIN32
		if (ppm_init == PPM_INIT_DUMP && ppm_recent.QuadPart < ppm_now.QuadPart)
#else
		if (ppm_init == PPM_INIT_DUMP && ppm_recent.tv_sec < ppm_now.tv_sec)
#endif
			return;
		ppm_recent = ppm_now;
		ppm_init = PPM_INIT_RUN;
		return;
	}
	nsamples += len / 2UL;
#ifdef _WIN32
	interval = ppm_now.QuadPart - ppm_recent.QuadPart;
	if (interval < (ppm_duration * perfFrequency.QuadPart))
		return;
#else
	interval = ppm_now.tv_sec - ppm_recent.tv_sec;
	if (interval < ppm_duration)
		return;
	interval *= 1000000000UL;
	interval += (int64_t)(ppm_now.tv_nsec - ppm_recent.tv_nsec);
#endif

	nsamples_total += nsamples;
	interval_total += interval;
	printf("real sample rate: %i current PPM: %i cumulative PPM: %i\n",
#ifdef _WIN32
		(int)get_rate(nsamples, interval, perfFrequency.QuadPart),
		ppm_report(nsamples, interval, perfFrequency.QuadPart),
		ppm_report(nsamples_total, interval_total, perfFrequency.QuadPart));

#else
		(int)get_rate(nsamples, interval, 1000000000UL),
		ppm_report(nsamples, interval, 1000000000UL),
		ppm_report(nsamples_total, interval_total, 1000000000UL));
#endif
	ppm_recent = ppm_now;
	nsamples = 0;
}

static void rtlsdr_callback(unsigned char *buf, uint32_t len, void *ctx)
{
	underrun_test(buf, len, 0);
	if (test_mode == PPM_BENCHMARK)
		ppm_test(len);
}

/* smallest band or band gap that tuner_benchmark() will notice */
static uint32_t max_step(uint32_t freq) {
	if (freq < 1000000U)
		return 10000U;
	if (freq > 100000000U)
		return 1000000U;
	return freq / 100U;
}

/* precision with which tuner_benchmark() will measure the edges of bands */
static uint32_t min_step(uint32_t freq) {
	return 100;
}

int confirm_pll_lock(uint32_t f)
{
	int i;
	for (i=0; i<20; i++) {
		if (rtlsdr_set_center_freq(dev, f) >= 0)
			return 1;
	}
	return 0;
}

/* returns last frequency before achieving status of 'lock' */
uint32_t coarse_search(uint32_t start, int lock)
{
	uint32_t f = start, f2;
	int status;
	while (f < SCAN_LIMIT) {
		if (do_exit)
			break;
		f2 = f + max_step(f);
		status = rtlsdr_set_center_freq(dev, f2) >= 0;
		if (!lock && !status)
			status = confirm_pll_lock(f2);
		if (status == lock)
			return f;
		f = f2;
	}
	return SCAN_LIMIT + 1;
}

/* returns frequency of a transition
 * must have one transition between start and start+step */
uint32_t fine_search(uint32_t start, uint32_t step)
{
	int low_status, mid_status, high_status;
	uint32_t f, stop;
	stop = start + step;
	f = start + step / 2;
	low_status = rtlsdr_set_center_freq(dev, start) >= 0;
	high_status = rtlsdr_set_center_freq(dev, stop) >= 0;
	if (low_status == high_status)
		return start;
	while (step > min_step(start)) {
		if (do_exit)
			break;
		mid_status = rtlsdr_set_center_freq(dev, f) >= 0;
		if (low_status == mid_status)
			start = f;
		else
			stop = f;
		step = stop - start;
		f = start + step / 2;
	}
	return f;
}

void tuner_benchmark(void)
{
	uint32_t f = 0, low_bound, high_bound;
	fprintf(stderr, "Testing tuner range. This may take a couple of minutes...\n");
	while (!do_exit) {
		/* find start of a band */
		f = coarse_search(f, 1);
		low_bound = fine_search(f, max_step(f));
		f += max_step(f);
		/* find stop of a band */
		f = coarse_search(f, 0);
		high_bound = fine_search(f, max_step(f));
		f += max_step(f);
		if (f > SCAN_LIMIT)
			break;
		fprintf(stderr, "Band: %u - %u Hz\n", low_bound, high_bound);
	}
}

int main(int argc, char **argv)
{
#ifndef _WIN32
	struct sigaction sigact;
#endif
	int n_read, r, opt, i;
	int sync_mode = 0;
	uint8_t *buffer;
	int dev_index = 0;
	int dev_given = 0;
	uint32_t out_block_size = DEFAULT_BUF_LENGTH;
	int count;
	int gains[100];

	while ((opt = getopt(argc, argv, "d:s:b:tp::Sh")) != -1) {
		switch (opt) {
		case 'd':
			dev_index = verbose_device_search(optarg);
			dev_given = 1;
			break;
		case 's':
			samp_rate = (uint32_t)atof(optarg);
			break;
		case 'b':
			out_block_size = (uint32_t)atof(optarg);
			break;
		case 't':
			test_mode = TUNER_BENCHMARK;
			break;
		case 'p':
			test_mode = PPM_BENCHMARK;
			if (optarg)
				ppm_duration = atoi(optarg);
			break;
		case 'S':
			sync_mode = 1;
			break;
		case 'h':
		default:
			usage();
			break;
		}
	}

	if(out_block_size < MINIMAL_BUF_LENGTH ||
	   out_block_size > MAXIMAL_BUF_LENGTH ){
		fprintf(stderr,
			"Output block size wrong value, falling back to default\n");
		fprintf(stderr,
			"Minimal length: %u\n", MINIMAL_BUF_LENGTH);
		fprintf(stderr,
			"Maximal length: %u\n", MAXIMAL_BUF_LENGTH);
		out_block_size = DEFAULT_BUF_LENGTH;
	}

	buffer = malloc(out_block_size * sizeof(uint8_t));

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
#else
	SetConsoleCtrlHandler( (PHANDLER_ROUTINE) sighandler, TRUE );
#endif
	count = rtlsdr_get_tuner_gains(dev, NULL);
	fprintf(stderr, "Supported gain values (%d): ", count);

	count = rtlsdr_get_tuner_gains(dev, gains);
	for (i = 0; i < count; i++)
		fprintf(stderr, "%.1f ", gains[i] / 10.0);
	fprintf(stderr, "\n");

	/* Set the sample rate */
	verbose_set_sample_rate(dev, samp_rate);

	if (test_mode == TUNER_BENCHMARK) {
		tuner_benchmark();
		goto exit;
	}

	/* Enable test mode */
	r = rtlsdr_set_testmode(dev, 1);

	/* Reset endpoint before we start reading from it (mandatory) */
	verbose_reset_buffer(dev);

	if ((test_mode == PPM_BENCHMARK) && !sync_mode) {
		fprintf(stderr, "Reporting PPM error measurement every %i seconds...\n", ppm_duration);
		fprintf(stderr, "Press ^C after a few minutes.\n");
	}

	if (test_mode == NO_BENCHMARK) {
		fprintf(stderr, "\nInfo: This tool will continuously"
				" read from the device, and report if\n"
				"samples get lost. If you observe no "
				"further output, everything is fine.\n\n");
	}

	if (sync_mode) {
		fprintf(stderr, "Reading samples in sync mode...\n");
		fprintf(stderr, "(Samples are being lost but not reported.)\n");
		while (!do_exit) {
			r = rtlsdr_read_sync(dev, buffer, out_block_size, &n_read);
			if (r < 0) {
				fprintf(stderr, "WARNING: sync read failed.\n");
				break;
			}

			if ((uint32_t)n_read < out_block_size) {
				fprintf(stderr, "Short read, samples lost, exiting!\n");
				break;
			}
			underrun_test(buffer, n_read, 1);
		}
	} else {
		fprintf(stderr, "Reading samples in async mode...\n");
		r = rtlsdr_read_async(dev, rtlsdr_callback, NULL,
				      0, out_block_size);
	}

	if (do_exit) {
		fprintf(stderr, "\nUser cancel, exiting...\n");
		fprintf(stderr, "Samples per million lost (minimum): %i\n", (int)(1000000L * dropped_samples / total_samples));
	}
	else
		fprintf(stderr, "\nLibrary error %d, exiting...\n", r);

exit:
	rtlsdr_close(dev);
	free (buffer);

	return r >= 0 ? r : -r;
}
