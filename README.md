# filter-c
Butterworth and chebyshev [lowpass, highpass, bandpass, bandstop] filters implementations in C.

based on algorithm from http://www.exstrom.com/journal/sigproc/

## Run example
```
make example
./example
```

## Steps to use a filter,
1. Create a filter object using `create_***_***_pass(params...)`
2. Use filter to filter incoming numbers one by one. The output is a double or float that can be specified in header.
3. After using the filter, release the filter using `free_***_***(filter)`.
