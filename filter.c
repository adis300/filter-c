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

double bw_low_pass(BWLowPass* filter, double x){
    int i;
    for(i=0; i<filter->n; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i] + x;
        x = filter->A[i]*(filter->w0[i] + 2.0*filter->w1[i] + filter->w2[i]);
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return x;
}
double bw_high_pass(BWHighPass* filter, double x){
    int i;
    for(i=0; i<filter->n; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i] + x;
        x = filter->A[i]*(filter->w0[i] - 2.0*filter->w1[i] + filter->w2[i]);
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return x;
}
double bw_band_pass(BWBandPass* filter, double x){
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
double bw_band_stop(BWBandStop* filter, double x){
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

CHELowPass* create_che_low_pass_filter(int n, double ep, double s, double f){
    CHELowPass* filter = (CHELowPass *) malloc(sizeof(CHELowPass));
    filter -> m = n/2;
    filter -> ep = 2.0/ep;
    filter -> A = (double *)malloc(filter -> m*sizeof(double));
    filter -> d1 = (double *)malloc(filter -> m*sizeof(double));
    filter -> d2 = (double *)malloc(filter -> m*sizeof(double));
    filter -> w0 = (double *)calloc(filter -> m, sizeof(double));
    filter -> w1 = (double *)calloc(filter -> m, sizeof(double));
    filter -> w2 = (double *)calloc(filter -> m, sizeof(double));

    double a = tan(M_PI*f/s);
    double a2 = a*a;
    double u = log((1.0+sqrt(1.0+ep*ep))/ep);
    double su = sinh(u/(double)n);
    double cu = cosh(u/(double)n);
    double b, c;
    
    int i;
    for(i=0; i < filter -> m; ++i){
        b = sin(M_PI*(2.0*i+1.0)/(2.0*n))*su;
        c = cos(M_PI*(2.0*i+1.0)/(2.0*n))*cu;
        c = b*b + c*c;
        s = a2*c + 2.0*a*b + 1.0;
        filter->A[i] = a2/(4.0*s); // 4.0
        filter->d1[i] = 2.0*(1-a2*c)/s;
        filter->d2[i] = -(a2*c - 2.0*a*b + 1.0)/s;
    }
    return filter;
}
CHEHighPass* create_che_high_pass_filter(int n, double ep, double s, double f){
    CHEHighPass* filter = (CHEHighPass *) malloc(sizeof(CHEHighPass));
    filter -> m = n/2;
    filter -> ep = 2.0/ep; // Used to normalize
    filter -> A = (double *)malloc(filter -> m*sizeof(double));
    filter -> d1 = (double *)malloc(filter -> m*sizeof(double));
    filter -> d2 = (double *)malloc(filter -> m*sizeof(double));
    filter -> w0 = (double *)calloc(filter -> m, sizeof(double));
    filter -> w1 = (double *)calloc(filter -> m, sizeof(double));
    filter -> w2 = (double *)calloc(filter -> m, sizeof(double));

    double a = tan(M_PI*f/s);
    double a2 = a*a;
    double u = log((1.0+sqrt(1.0+ep*ep))/ep);
    double su = sinh(u/(double)n);
    double cu = cosh(u/(double)n);
    double b, c;
    
    int i;
    for(i=0; i < filter -> m; ++i){
        b = sin(M_PI*(2.0*i+1.0)/(2.0*n))*su;
        c = cos(M_PI*(2.0*i+1.0)/(2.0*n))*cu;
        c = b*b + c*c;
        s = a2 + 2.0*a*b + c;
        filter->A[i] = 1.0/(4.0*s); // 4.0
        filter->d1[i] = 2.0*(c-a2)/s;
        filter->d2[i] = -(a2 - 2.0*a*b + c)/s;
    }
    return filter;
}
CHEBandPass* create_che_band_pass_filter(int n, double ep, double s, double f1, double f2){

}
CHEBandStop* create_che_band_stop_filter(int n, double ep, double s, double f1, double f2){

}

void free_che_low_pass(CHELowPass* filter){
    free(filter -> A);
    free(filter -> d1);
    free(filter -> d2);
    free(filter -> w0);
    free(filter -> w1);
    free(filter -> w2);
    free(filter);
}
void free_che_high_pass(CHEHighPass* filter){
    free(filter -> A);
    free(filter -> d1);
    free(filter -> d2);
    free(filter -> w0);
    free(filter -> w1);
    free(filter -> w2);
    free(filter);
}
void free_che_band_pass(CHEBandPass* filter){
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
void free_che_band_stop(CHEBandStop* filter){
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

double che_low_pass(CHELowPass* filter, double x){
    int i;
    for(i=0; i<filter->m; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i] + x;
        x = filter->A[i]*(filter->w0[i] + 2.0*filter->w1[i] + filter->w2[i]);
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return filter->ep*x;
}
double che_high_pass(CHEHighPass* filter, double x){
    int i;
    for(i=0; i<filter->m; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i] + x;
        x = filter->A[i]*(filter->w0[i] - 2.0*filter->w1[i] + filter->w2[i]);
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return filter->ep*x;
}
double che_band_pass(CHEBandPass* filter, double input){

}
double che_band_stop(CHEBandStop* filter, double input){
    
}