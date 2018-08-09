/*
 * filter.c
 *
 * Copyright (c) 2018 Disi A
 * 
 * Author: Disi A
 * Email: adis@live.cn
 *  https://www.mathworks.com/matlabcentral/profile/authors/3734620-disi-a
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "filter.h"

BWLowPass* create_bw_low_pass_filter(int order, double s, double f) {
    BWLowPass* filter = (BWLowPass *) malloc(sizeof(BWLowPass));
    filter -> n = order/2;
    filter -> A = (double *)malloc(filter -> n*sizeof(double));
    filter -> d1 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d2 = (double *)malloc(filter -> n*sizeof(double));
    filter -> w0 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w1 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w2 = (double *)calloc(filter -> n, sizeof(double));

    double a = tan(M_PI * f/s);
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
BWHighPass* create_bw_high_pass_filter(int order, double s, double f){
    BWHighPass* filter = (BWHighPass *) malloc(sizeof(BWHighPass));
    filter -> n = order/2;
    filter -> A = (double *)malloc(filter -> n*sizeof(double));
    filter -> d1 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d2 = (double *)malloc(filter -> n*sizeof(double));
    filter -> w0 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w1 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w2 = (double *)calloc(filter -> n, sizeof(double));

    double a = tan(M_PI * f/s);
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
BWBandPass* create_bw_band_pass_filter(int order, double s, double f1, double f2){
    BWBandPass* filter = (BWBandPass *) malloc(sizeof(BWBandPass));
    filter -> n = order/4;
    filter -> A = (double *)malloc(filter -> n*sizeof(double));
    filter -> d1 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d2 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d3 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d4 = (double *)malloc(filter -> n*sizeof(double));

    filter -> w0 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w1 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w2 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w3 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w4 = (double *)calloc(filter -> n, sizeof(double));

  
    double a = cos(M_PI*(f1+f2)/s)/cos(M_PI*(f1-f2)/s);
    double a2 = a*a;
    double b = tan(M_PI*(f1-f2)/s);
    double b2 = b*b;
    double r;
    int i;
    for(i=0; i<filter->n; ++i){
        r = sin(M_PI*(2.0*i+1.0)/(4.0*filter->n));
        s = b2 + 2.0*b*r + 1.0;
        filter->A[i] = b2/s;
        filter->d1[i] = 4.0*a*(1.0+b*r)/s;
        filter->d2[i] = 2.0*(b2-2.0*a2-1.0)/s;
        filter->d3[i] = 4.0*a*(1.0-b*r)/s;
        filter->d4[i] = -(b2 - 2.0*b*r + 1.0)/s;
    }
    return filter;
}
BWBandStop* create_bw_band_stop_filter(int order, double s, double f1, double f2){
    BWBandStop* filter = (BWBandStop *) malloc(sizeof(BWBandStop));
    filter -> n = order/4;
    filter -> A = (double *)malloc(filter -> n*sizeof(double));
    filter -> d1 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d2 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d3 = (double *)malloc(filter -> n*sizeof(double));
    filter -> d4 = (double *)malloc(filter -> n*sizeof(double));

    filter -> w0 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w1 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w2 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w3 = (double *)calloc(filter -> n, sizeof(double));
    filter -> w4 = (double *)calloc(filter -> n, sizeof(double));

    double a = cos(M_PI*(f1+f2)/s)/cos(M_PI*(f1-f2)/s);
    double a2 = a*a;
    double b = tan(M_PI*(f1-f2)/s);
    double b2 = b*b;
    double r;

    int i;
    for(i=0; i<filter->n; ++i){
        r = sin(M_PI*(2.0*i+1.0)/(4.0*filter->n));
        s = b2 + 2.0*b*r + 1.0;
        filter->A[i] = 1.0/s;
        filter->d1[i] = 4.0*a*(1.0+b*r)/s;
        filter->d2[i] = 2.0*(b2-2.0*a2-1.0)/s;
        filter->d3[i] = 4.0*a*(1.0-b*r)/s;
        filter->d4[i] = -(b2 - 2.0*b*r + 1.0)/s;
    }
    filter->r = 4.0*a;
    filter->s = 4.0*a2+2.0;
    return filter;
}

void free_bw_low_pass(BWLowPass* filter){
    free(filter -> A);
    free(filter -> d1);
    free(filter -> d2);
    free(filter -> w0);
    free(filter -> w1);
    free(filter -> w2);
    free(filter);
}
void free_bw_high_pass(BWHighPass* filter){
    free(filter -> A);
    free(filter -> d1);
    free(filter -> d2);
    free(filter -> w0);
    free(filter -> w1);
    free(filter -> w2);
    free(filter);
}
void free_bw_band_pass(BWBandPass* filter){
    free(filter -> A);
    free(filter -> d1);
    free(filter -> d2);
    free(filter -> d3);
    free(filter -> d4);
    free(filter -> w0);
    free(filter -> w1);
    free(filter -> w2);
    free(filter -> w3);
    free(filter -> w4);
    free(filter);
}
void free_bw_band_stop(BWBandStop* filter){
    free(filter -> A);
    free(filter -> d1);
    free(filter -> d2);
    free(filter -> d3);
    free(filter -> d4);
    free(filter -> w0);
    free(filter -> w1);
    free(filter -> w2);
    free(filter -> w3);
    free(filter -> w4);
    free(filter);
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
double band_pass(BWBandPass* filter, double x){
    int i;
    for(i=0; i<filter->n; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i]+ filter->d3[i]*filter->w3[i]+ filter->d4[i]*filter->w4[i] + x;
        x = filter->A[i]*(filter->w0[i] - 2.0*filter->w2[i] + filter->w4[i]);
        filter->w4[i] = filter->w3[i];
        filter->w3[i] = filter->w2[i];
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return x;
}
double band_stop(BWBandStop* filter, double x){
    int i;
    for(i=0; i<filter->n; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i]+ filter->d3[i]*filter->w3[i]+ filter->d4[i]*filter->w4[i] + x;
        x = filter->A[i]*(filter->w0[i] - filter->r*filter->w1[i] + filter->s*filter->w2[i]- filter->r*filter->w3[i] + filter->w4[i]);
        filter->w4[i] = filter->w3[i];
        filter->w3[i] = filter->w2[i];
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return x;
}