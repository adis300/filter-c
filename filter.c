#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "filter.h"

BWLowPass* create_bw_low_pass_filter(int order, double sampling_frequency, int half_power_frequency) {
    BWLowPass* filter = (BWLowPass *) malloc(sizeof(BWLowPass));
    filter -> n = order/2;
    filter -> A = (double *)malloc(filter -> n*sizeof(double));
    filter -> d1 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d2 = (double *)malloc(filter -> n*sizeof(double));
    filter -> w0 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w1 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w2 = (double *)calloc(filter -> n, sizeof(double));

    double s = sampling_frequency;
    double a = tan(M_PI * half_power_frequency/s);
    double a2 = a * a;
    double r;
    
    int i;
    
    for(i=0; i < filter -> n; ++i){
        r = sin(M_PI*(2.0*i+1.0)/(4.0*filter -> n));
        s = a2 + 2.0*a*r + 1.0;
        filter -> A[i] = a2/s;
        filter -> d1[i] = 2.0*(1-a2)/s;
        filter -> d2[i] = -(a2 - 2.0*a*r + 1.0)/s;
    }
    return filter;
}
BWHighPass* create_bw_high_pass_filter(int order, double sampling_frequency, int half_power_frequency){
    BWHighPass* filter = (BWHighPass *) malloc(sizeof(BWHighPass));
    filter -> n = order/2;
    filter -> A = (double *)malloc(filter -> n*sizeof(double));
    filter -> d1 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d2 = (double *)malloc(filter -> n*sizeof(double));
    filter -> w0 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w1 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w2 = (double *)calloc(filter -> n, sizeof(double));

    double s = sampling_frequency;
    double a = tan(M_PI * half_power_frequency/s);
    double a2 = a * a;
    double r;
    
    int i;

    for(i=0; i < filter -> n; ++i){
        r = sin(M_PI*(2.0*i+1.0)/(4.0*filter -> n));
        s = a2 + 2.0*a*r + 1.0;
        filter -> A[i] = 1.0/s;
        filter -> d1[i] = 2.0*(1-a2)/s;
        filter -> d2[i] = -(a2 - 2.0*a*r + 1.0)/s;
    }
    return filter;
}
BWBandPass* create_bw_band_pass_filter(int order, double sampling_frequency, int half_power_frequency1, int half_power_frequency2){
    return NULL;
}
BWBandStop* create_bw_band_stop_filter(int order, double sampling_frequency, int half_power_frequency1, int half_power_frequency2){
    return NULL;
}

void free_bw_lphp(BWLowPass* filter){
    free(filter -> d1);
    free(filter -> d2);
    free(filter -> w0);
    free(filter -> w1);
    free(filter -> w2);
    free(filter);
}

void free_bw_low_pass(BWLowPass* filter){
    free_bw_lphp(filter);
}
void free_bw_high_pass(BWHighPass* filter){
    free_bw_lphp(filter);
}
void free_bw_band_pass(BWBandPass* filter){

}
void free_bw_band_stop(BWBandStop* filter){

}

double low_pass(BWLowPass* filter, double x){
    int i;
    for(i=0; i<filter->n; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i] + x;
        x = filter->A[i]*(filter->w0[i] + 2.0*filter->w1[i] + filter->w2[i]);
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return x;
}
double high_pass(BWHighPass* filter, double x){
    int i;
    for(i=0; i<filter->n; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i] + x;
        x = filter->A[i]*(filter->w0[i] - 2.0*filter->w1[i] + filter->w2[i]);
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return x;
}
double band_pass(BWBandPass* filter, double input){

}
double band_stop(BWBandStop* filter, double input){

}