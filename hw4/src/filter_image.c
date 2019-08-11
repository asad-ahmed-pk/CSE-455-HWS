// filter_image.c
// Image filtering hw 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"

#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    float sum = 0.0;
    int x, y, c;

    // calculate sum
    for (x = 0; x < im.w; x++)
    {
        for (y = 0; y < im.h; y++)
        {
            for (c = 0; c < im.c; c++)
            {
                sum += (get_pixel(im, x, y, c));
            }
        }
    }

    // to avoid divide by zero
    if (sum == 0) {
        return;
    }

    // divide each pixel value by the sum
    float p = 0.0;
    for (x = 0; x < im.w; x++)
    {
        for (y = 0; y < im.h; y++)
        {
            for (c = 0; c < im.c; c++)
            {
                p = get_pixel(im, x, y, c);
                set_pixel(im, x, y, c, p / sum);
            }
        }
    }
}

image make_box_filter(int w)
{
    // TODO
    image filter = make_image(w, w, 1);
    
    // set all to 1 / wxw
    for (int x = 0; x < filter.w; x++) {
        for (int y = 0; y < filter.h; y++) {
            set_pixel(filter, x, y, 0, 1.0 / (w*w));
        }
    }

    return filter;
}

// Get the sum for the given pixel point (ix, iy) in the image im for the given filter for the given channel
// ix, iy, ic: the image pixel coordinate and channel
// fc: the filter channel
float get_sum_for_convolution_at(image im, int ix, int iy, int ic, image filter, int fc)
{
    // calculate the start and end indices (bounding box for the neighbourhood)
    int x_start = ix - (int)(filter.w / 2) ;
    int x_end = ix + (int)(filter.w / 2);
    
    int y_start = iy - (int)(filter.h / 2);
    int y_end = iy + (int)(filter.h / 2);

    // image indices
    int x, y;

    // filter indices
    int fx = 0; int fy = 0;       

    float sum = 0.0;
    float value = 0.0;

    for (y = y_start; y <= y_end; y++) 
    {
        for (x = x_start; x <= x_end; x++)  
        {
            if (x < 0 || x >= im.w || y < 0 || y >= im.h) {
                value = 0.0;
            }
            else {
                value = (get_pixel(im, x, y, ic) * get_pixel(filter, fx, fy, fc));
            }

            sum += value;

            fx = (fx + 1) % filter.w;
        }

        fy = (fy + 1) % filter.h;
    }

    return sum;
}

// Convolve the image at
void convole_at(image im, image result, image filter, int ix, int iy, int perserve)
{
    float sum = 0.0;

    // same channels in image and filter
    if (im.c == filter.c)
    {
        for (int c = 0; c < im.c; c++) 
        {
            sum = get_sum_for_convolution_at(im, ix, iy, c, filter, c);

            if (perserve) {
                set_pixel(result, ix, iy, c, sum);
            }
        }

        if (!perserve) {
            set_pixel(result, ix, iy, 0, sum);
        }
    }
    else
    {
        for (int c = 0; c < im.c; c++) 
        {
            sum = get_sum_for_convolution_at(im, ix, iy, c, filter, 0);

            if (perserve) {
                set_pixel(result, ix, iy, c, sum);
            }
        }

        if (!perserve) {
            set_pixel(result, ix, iy, 0, sum);
        }
    }
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    
    // filter must have 1 channel or have same channels as image
    assert(filter.c == 1 || filter.c == im.c);

    // create the result image
    image result = make_image(im.w, im.h, preserve ? im.c : 1);

    float h = 0.0;
    float g = 0.0;
    float v = 0.0;

    for (int i = 0; i < im.w; i++)
    {
        for (int j = 0; j < im.h; j++)
        {
            for (int c = 0; c < im.c; c++)
            {
                // if preserving - need to reset to 0 so sum is calculated for this channel
                if (preserve) {
                    v = 0.0;
                }

                // m, n - the filter indices
                for (int m = 0; m < filter.w; m++)
                {
                    for (int n = 0; n < filter.h; n++)
                    {
                        h = get_pixel(im, i - m, j - n, c);
                        g = get_pixel(filter, m, n, preserve ? (filter.c > 1 ? c : 0) : 0);
                        v += h * g;
                    }
                }

                set_pixel(result, i, j, preserve ? c : 0, v);
            }

            v = 0.0;
        }
    }

    return result;
}

image make_highpass_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);
    float matrix[9] = {0, -1, 0, -1, 4, -1, 0, -1, 0};

    memcpy(filter.data, matrix, sizeof(matrix));

    return filter;
}

image make_sharpen_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);
    float matrix[9] = {0, -1, 0, -1, 5, -1, 0, -1, 0};

    memcpy(filter.data, matrix, sizeof(matrix));

    return filter;
}

image make_emboss_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);
    float matrix[9] = {-2, -1, 0, -1, 1, 1, 0, 1, 2};

    memcpy(filter.data, matrix, sizeof(matrix));

    return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: Box filter should use the preserve option because the requirement is to smooth the image. This operation needs to be applied
// to all 3 channels to get a valid output.
// Similarly, the output is required to be even for the emboss and sharpen filter - it needs to be applied to all 3 channels.
// The highpass filter only allows for high frequencies to be included in the final image - hence

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Yes, this is needed for the high-pass filter as it leads to an overflow of the pixel values and needs to be clamped.

image make_gaussian_filter(float sigma)
{
    // TODO
    int w = sigma == 0 ? 1 : 6 * sigma;
    if (w % 2 == 0) {
        w = w + 1;
    }

    image kernel = make_image(w, w, 1);

    int x, y;
    int offset = (int)(w / 2);
    float G = 0;

    for (y = -offset; y < w - offset; y++)
    {
        for (x = -offset; x < w - offset; x++)
        {
            G = (1.0 / (TWOPI * sigma * sigma)) * (expf(-(x * x + y * y) / (2 * sigma * sigma)));
            set_pixel(kernel, x + offset, y + offset, 0, G);
        }
    }

    l1_normalize(kernel);

    return kernel;
}

image add_image(image a, image b)
{
    // TODO
    
    // assert same dimensions
    assert(a.w == b.w && a.h == b.h && a.c == b.c);

    image result = make_image(a.w, a.h, a.c);

    float x, y, c;
    for (y = 0; y < a.h; y++)
    {
        for (x = 0; x < a.w; x++)
        {
            for (c = 0; c < a.c; c++)
            {
                set_pixel(result, x, y, c, get_pixel(a, x, y, c) + get_pixel(b, x, y, c));
            }
        }
    }

    return result;
}

image sub_image(image a, image b)
{
    // TODO

    image result = make_image(a.w, a.h, a.c);

    float x, y, c;
    for (y = 0; y < a.h; y++)
    {
        for (x = 0; x < a.w; x++)
        {
            for (c = 0; c < a.c; c++)
            {
                set_pixel(result, x, y, c, get_pixel(a, x, y, c) - get_pixel(b, x, y, c));
            }
        }
    }

    return result;
}

image make_gx_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);
    float matrix[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};

    memcpy(filter.data, matrix, sizeof(matrix));

    return filter;
}

image make_gy_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);
    float matrix[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

    memcpy(filter.data, matrix, sizeof(matrix));

    return filter;
}

void feature_normalize(image im)
{    
    // TODO

    // finding the min and max
    float min = im.data[0];
    float max = min;

    for (int i = 0; i < im.c * im.w * im.h; i++) 
    {
        if (min > im.data[i]) min = im.data[i];
        if (max < im.data[i]) max = im.data[i];
    }

    float range = max - min;

    // safeguard against divide by zero
    if (range == 0) 
    {
        for (int i = 0; i < im.c * im.w * im.h; i++) {
            im.data[i] = 0;
        }
    }

    for (int i = 0; i < im.c * im.w * im.h; i++)  {
        im.data[i] = (im.data[i] - min) / range;
    }
}

image* sobel_image(image im)
{
    // TODO
    image* images = calloc(2, sizeof(image));

    // calculate Gx, and Gy
    image gx_filter = make_gx_filter();
    image Gx = convolve_image(im, gx_filter, 0);
    
    image gy_filter = make_gx_filter();
    image Gy = convolve_image(im, gy_filter, 0);

    // calculate gradient
    image gradient = make_image(im.w, im.h, 1);
    for (int i = 0; i < im.w * im.h; i++) {
        gradient.data[i] = sqrtf(Gx.data[i] * Gx.data[i] + Gy.data[i] * Gy.data[i]);
    }

    // calculate theta (direction)
    image theta = make_image(im.w, im.h, 1);
    for (int i = 0; i < im.w * im.h; i++) {
        theta.data[i] = atan2f(Gy.data[i], Gx.data[i]);
    }

    images[0] = gradient;
    images[1] = theta;

    free_image(gx_filter);
    free_image(gy_filter);
    free_image(Gx);
    free_image(Gy);

    return images;
}

image colorize_sobel(image im)
{
    // TODO
    image* sobel_result = sobel_image(im);
    image magnitude = sobel_result[0];
    image direction = sobel_result[1];

    // need to normalize the results
    feature_normalize(magnitude);
    feature_normalize(direction);

    // create hsv image
    image hsv = make_image(im.w, im.h, 3);

    // angle for the hue
    float x, y;
    for (y = 0; y < im.h; y++)
    {
        for (x = 0; x < im.w; x++)
        {
            // direction for hue
            set_pixel(hsv, x, y, 0, get_pixel(direction, x, y, 0));

            // magnitude for saturation and value
            set_pixel(hsv, x, y, 1, get_pixel(magnitude, x, y, 0));
            set_pixel(hsv, x, y, 2, get_pixel(magnitude, x, y, 0));
        }
    }

    hsv_to_rgb(hsv);

    return hsv;
}
