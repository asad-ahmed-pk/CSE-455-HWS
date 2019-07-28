#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: optional, make separable 1d Gaussian.
    return make_image(1,1,1);
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(1){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.
        return copy_image(im);
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    // TODO: calculate structure matrix for im.
    image S = make_image(im.w, im.h, 3);

    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();

    image ix = convolve_image(im, gx_filter, 0);
    image iy = convolve_image(im, gy_filter, 0);

    image I = make_image(im.w, im.h, 3);

    int x, y;
    for (y = 0; y < im.h; y++)
    {
        for (x = 0; x < im.w; x++)
        {
            // set the values if IxIx, IyIy, and IxIy in 3 separate channels
            set_pixel(I, x, y, 0, get_pixel(ix, x, y, 0) * get_pixel(ix, x, y, 0));
            set_pixel(I, x, y, 1, get_pixel(iy, x, y, 0) * get_pixel(iy, x, y, 0));
            set_pixel(I, x, y, 2, get_pixel(ix, x, y, 0) * get_pixel(iy, x, y, 0));            
        }
    }

    // S is the smoothed image of I (weighted sum)
    S = smooth_image(I, sigma);

    free_image(I);
    free_image(ix);
    free_image(iy);
    free_image(gx_filter);
    free_image(gy_filter);

    save_image(S, "structure_matrix");

    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);

    static float alpha = 0.06;

    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.

    int x, y;
    float value = 0.0;
    float det = 0.0;
    float trace = 0.0;
    float ixx, ixy, iyy;

    for (y = 0; y < S.h; y++)
    {
        for (x = 0; x < S.w; x++)
        {
            ixx = get_pixel(S, x, y, 0);
            iyy = get_pixel(S, x, y, 1);
            ixy = get_pixel(S, x, y, 2);

            det = ixx * iyy - ixy * ixy;
            trace = ixx + iyy;

            value = det - alpha * (trace * trace);

            set_pixel(R, x, y, 0, value);
        }
    }

    save_image(R, "cornerness_response");

    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    static float LOW_RESPONSE = -999999;

    int x, y;
    int x_inner, y_inner;
    int x_min, x_max, y_min, y_max;

    for (y = 0; y < im.h; y++)
    {
        for (x = 0; x < im.w; x++)
        {
            x_min = MAX(0, x - w);
            x_max = MIN(im.w - 1, x + w);
            y_min = MAX(0, y - w);
            y_max = MIN(im.h - 1, y + w);

            for (x_inner = x_min; x_inner <= x_max; x_inner++)
            {
                for (y_inner = y_min; y_inner <= y_max; y_inner++)
                {
                    if (get_pixel(im, x_inner, y_inner, 0) > get_pixel(im, x, y, 0)) {
                        set_pixel(r, x, y, 0, LOW_RESPONSE);
                    }
                }
            }
        }
    }

    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    save_image(Rnms, "RNMS");

    //TODO: count number of responses over threshold
    int count = 0;

    for (int y = 0; y < Rnms.h; y++)
    {
        for (int x = 0; x < Rnms.w; x++)
        {
            if (get_pixel(Rnms, x, y, 0) > thresh) {
                count++;
            }
        }
    }
    
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));

    //TODO: fill in array *d with descriptors of corners, use describe_index.
    int index = 0;
    for (int i = 0; i < Rnms.w * Rnms.h; i++) 
    {
        if (Rnms.data[i] > thresh) {
            d[index] = describe_index(im, i);
            index++;
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);

    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
