// resize_image.c
// Image resizing Homework 1

#include <math.h>
#include <stdio.h>

#include "image.h"

// min and max functions
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// Interpolate the given coordinates (x, y) from the image im and return the pixel value
float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    int x1 = round(x + 0.5);
    int y1 = round(y + 0.5);

    return get_pixel(im, x1, y1, c);
}

image nn_resize(image im, int w, int h)
{
    // TODO
    image img_new = make_image(w, h, im.c);

    float width_ratio = (float)im.w / (float)w;
    float height_ratio = (float)im.h / (float)h;
    float value = 0.0;

    int x1, y1;

    // resize the image
    for (int y = 0; y < img_new.h; y++)
    {
        for (int x = 0; x < img_new.w; x++)
        {
            for (int c = 0; c < img_new.c; c++)
            {
                x1 = x * width_ratio;
                y1 = y * height_ratio;
                value = nn_interpolate(im, x1, y1, c);

                set_pixel(img_new, x, y, c, value);
            }
        }
    }

    return img_new;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    // find the 4 pixels in the image that form the bounding box around q (x,y)
    float v1_x = floorf(x);
    float v1_y = floorf(y);

    float v2_x = ceilf(x);
    float v2_y = floorf(y);

    float v3_x = floorf(x);
    float v3_y = ceilf(y);

    float v4_x = ceilf(x);
    float v4_y = ceilf(y);

    // get their pixel values
    float v1 = get_pixel(im, floor(v1_x), floor(v1_y), c);
    float v2 = get_pixel(im, ceil(v2_x), floor(v2_y), c);
    float v3 = get_pixel(im, floor(v3_x), ceil(v3_y), c);
    float v4 = get_pixel(im, ceil(v4_x), ceil(v4_y), c);

    // calculate distances
    float d1 = x - v1_x;
    float d2 = v2_x - x;
    float d3 = y - v2_y;
    float d4 = v4_y - y;

    // calculate qs
    float q1 = v1 * d2 + v2 * d1;
    float q2 = v3 * d2 + v4 * d1;
    float q = q1 * d4 + q2 * d3;

    return q;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    image image_new = make_image(w, h, im.c);

    float width_ratio = (float)im.w / (float)w;
    float height_ratio = (float)im.h / (float)h;
    float value = 0.0;

    float qx, qy;

    // resize the image
    for (int y = 0; y < image_new.h; y++)
    {
        for (int x = 0; x < image_new.w; x++)
        {
            for (int c = 0; c < image_new.c; c++)
            {
                qx = (x + 0.5) * width_ratio;
                qy = (y + 0.5) * height_ratio;
                value = bilinear_interpolate(im, qx, qy, c);

                set_pixel(image_new, x, y, c, value);
            }
        }
    }

    return image_new;
}
