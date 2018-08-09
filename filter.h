/*
 * filter.h
 *
 * Copyright (c) 2018 Disi A
 * 
 * Author: Disi A
 * Email: adis@live.cn
 *  https://www.mathworks.com/matlabcentral/profile/authors/3734620-disi-a
 */
#ifndef filter_h
#define filter_h

#if __cplusplus
extern "C"{
#endif

#define PRECISION double

typedef struct {
    int n;
	double *A;
    double *d1;
    double *d2;
    double *w0;
    double *w1;
    double *w2;
} BWLowPass;
typedef BWLowPass BWHighPass;

typedef struct {
    int n;
	double *A;
    double *d1;
    double *d2;
    double *d3;
    double *d4;
    double *w0;
    double *w1;
    double *w2;
    double *w3;
    double *w4;
} BWBandPass;

typedef struct {
    int n;
    double r;
    double s;
	double *A;
    double *d1;
    double *d2;
    double *d3;
    double *d4;
    double *w0;
    double *w1;
    double *w2;
    double *w3;
    double *w4;
} BWBandStop;

typedef struct {
    int n;
    double ep;
	double *A;
    double *d1;
    double *d2;
    double *w0;
    double *w1;
    double *w2;
} CHELowPass;
typedef CHELowPass CHEHighPass;

typedef struct {
    int n;
    double ep;
	double *A;
    double *d1;
    double *d2;
    double *d3;
    double *d4;
    double *w0;
    double *w1;
    double *w2;
    double *w3;
    double *w4;
} CHEBandPass;

typedef struct {
    int n;
    double ep;
    double r;
    double s;
	double *A;
    double *d1;
    double *d2;
    double *d3;
    double *d4;
    double *w0;
    double *w1;
    double *w2;
    double *w3;
    double *w4;
} CHEBandStop;

BWLowPass* create_bw_low_pass_filter(int order, double sampling_frequency, double half_power_frequency);
BWHighPass* create_bw_high_pass_filter(int order, double sampling_frequency, double half_power_frequency);
BWBandPass* create_bw_band_pass_filter(int order, double sampling_frequency, double half_power_frequency1, double half_power_frequency2);
BWBandStop* create_bw_band_stop_filter(int order, double sampling_frequency, double half_power_frequency1, double half_power_frequency2);

void free_bw_low_pass(BWLowPass* filter);
void free_bw_high_pass(BWHighPass* filter);
void free_bw_band_pass(BWBandPass* filter);
void free_bw_band_stop(BWBandStop* filter);

double low_pass(BWLowPass* filter, double input);
double high_pass(BWHighPass* filter, double input);
double band_pass(BWBandPass* filter, double input);
double band_stop(BWBandStop* filter, double input);

#if __cplusplus
}
#endif

#endif /* filter_h */