# testing python script

from uwimg import *

# 1. Image resizing
im = load_image("data/dogsmall.jpg")
a = nn_resize(im, im.w*4, im.h*4)
save_image(a, "dog4x-nn")

b = bilinear_resize(im, im.w*4, im.h*4)
save_image(b, "dog4x-bl")

c = nn_resize(im, im.w//7, im.h//7)
save_image(a, "dog7th-bl")

# 2. Convolutions

im = load_image("data/dog.jpg")

f = make_box_filter(7)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-box7")
thumb = nn_resize(blur, blur.w//7, blur.h//7)
save_image(thumb, "dogthumb")

f2 = make_highpass_filter()
blur = convolve_image(im, f2, 0)
clamp_image(blur)
save_image(blur, "dog-highpass")

f3 = make_sharpen_filter()
blur = convolve_image(im, f3, 1)
save_image(blur, "dog-sharpen")
clamp_image(blur)

f4 = make_emboss_filter()
blur = convolve_image(im, f4, 1)
save_image(blur, "dog-emboss")
clamp_image(blur)

# 2.3 Gaussian Kernel

im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-gauss2")

# 3. Hybrid Images
im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
lfreq = convolve_image(im, f, 1)
hfreq = im - lfreq
reconstruct = lfreq + hfreq
save_image(lfreq, "low-frequency")
save_image(hfreq, "high-frequency")
save_image(reconstruct, "reconstruct")

# creating Ronbledore
ron = load_image("data/ron.png")
dumbledore = load_image("data/dumbledore.png")

lfreq_ron = convolve_image(ron, f, 1)
hfreq_dumbledore = dumbledore - convolve_image(dumbledore, f, 1)
ronbledore = lfreq_ron + hfreq_dumbledore
save_image(ronbledore, "ronbledore")

# 4. Sobel Filters
im = load_image("data/dog.jpg")
res = sobel_image(im)
mag = res[0]
feature_normalize(mag)
save_image(mag, "magnitude")

col_sobel = colorize_sobel(im)
feature_normalize(col_sobel)
save_image(col_sobel, "sobel-color")
