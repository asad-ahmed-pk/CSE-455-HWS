#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "image.h"

// private functions
int get_index(image im, int x, int y, int c, unsigned short should_clamp);

float get_pixel(image im, int x, int y, int c) {
    return im.data[get_index(im, x, y, c, 1)];
}

void set_pixel(image im, int x, int y, int c, float v) {
    int index = get_index(im, x, y, c, 0);
    if (index < 0) return;
    im.data[index] = v;
}

image copy_image(image im) {
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, im.h * im.w * im.c);
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);

    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++) 
        {
            float r = get_pixel(im, x, y, 0);
            float g = get_pixel(im, x, y, 1);
            float b = get_pixel(im, x, y, 2);
            float v = 0.299 * r + 0.587 * g + 0.114 * b;

            set_pixel(gray, x, y, 0, v);
        }
    }

    return gray;
}

void shift_image(image im, int c, float v)
{
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++) {
            float orig_value = get_pixel(im, x, y, c);
            set_pixel(im, x, y, c, orig_value + v);
        }
    }
}

// extra credit
void scale_image(image im, int c, float v)
{
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++) {
            float orig_value = get_pixel(im, x, y, c);
            set_pixel(im, x, y, c, orig_value * v);
        }
    }
}

void clamp_image(image im)
{
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++)
        {
            for (int c = 0; c < im.c; c++)
            {
                float value = get_pixel(im, x, y, c);
                if (value < 0.0) value = 0.0;
                else if (value >= 1.0) value = 1.0;
                set_pixel(im, x, y, c, value);
            }
        }
    }
}

// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++)
        {
            float V = three_way_max(get_pixel(im, x, y, 0), get_pixel(im, x, y, 1), get_pixel(im, x, y, 2));
            float m = three_way_min(get_pixel(im, x, y, 0), get_pixel(im, x, y, 1), get_pixel(im, x, y, 2));
            float C = V - m;

            float S = (V <= 0.0 ? 0.0 : C / V);

            float h_prime = 0.0;
            if (C > 0.0) 
            {
                if (V == get_pixel(im, x, y, 0)) {
                    h_prime = (get_pixel(im, x, y, 1) - get_pixel(im, x, y, 2)) / C;
                }
                else if (V == get_pixel(im, x, y, 1)) {
                    h_prime = ((get_pixel(im, x, y, 2) - get_pixel(im, x, y, 0)) / C) + 2;
                }
                else {
                    h_prime = ((get_pixel(im, x, y, 0) - get_pixel(im, x, y, 1)) / C) + 4;
                }
            }

            float H = (h_prime <= 0.0) ? (h_prime / 6.0) + 1 : h_prime / 6.0;

            set_pixel(im, x, y, 0, H);
            set_pixel(im, x, y, 1, S);
            set_pixel(im, x, y, 2, V);
        }
    }
}

void hsv_to_rgb(image im)
{
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++)
        {
            int H = (int)(get_pixel(im, x, y, 0) * 360.0);
            float S = get_pixel(im, x, y, 1);
            float V = get_pixel(im, x, y, 2);

            float C = V * S;
            float X = C * (1 - abs(((H / 60) % 2) - 1));
            float m = V - C;

            float R = 0.0;
            float G = 0.0;
            float B = 0.0;

            if (H >= 0 && H < 60) {
                R = C; G = X;
            }
            else if (H >= 60 && H < 120) {
                R = X; G = C;
            }
            else if (H >= 120 && H < 180) {
                G = C; B = X;
            }
            else if (H >= 180 && H < 240) {
                G = X; B = C;
            }
            else if (H >= 240 && H < 300) {
                R = X; B = C;
            }
            else {
                R = C; B = X;
            }

            R += m; G += m; B += m;

            set_pixel(im, x, y, 0, R);
            set_pixel(im, x, y, 1, G);
            set_pixel(im, x, y, 2, B);
        }
    }
}

// Returns the index for the given x, y, c values (row, col, channel)
// If should_clamp is set to 1, index values are clamped. Otherwise returns -1 on out of bounds.
int get_index(image im, int x, int y, int c, unsigned short should_clamp) 
{
    // bounds checking for x (col)
    if (x >= im.w) 
    {
        if (should_clamp) {
            return ((im.w - 1) + y * (im.w) + c * im.w * im.h);
        }
        else {
            return -1;
        }
    }
    else if (x < 0)
    {
        if (should_clamp) {
            return (0 + y * (im.w) + c * im.w * im.h);
        }
        else {
            return -1;
        }
    }

    //bounds checking for y (row)
    if (y >= im.h)
    {
        if (should_clamp) {
            return (x + ((im.h - 1) * (im.w)) + c * im.w * im.h);
        }
        else {
            return -1;
        }
    }
    else if (y < 0) 
    {
        if (should_clamp) {
            return (x + 0 * (im.w) + c * im.w * im.h);
        }
        else {
            return -1;
        }
    }

    // valid index
    return (x + (y * (im.w)) + (c * im.w * im.h));
}

