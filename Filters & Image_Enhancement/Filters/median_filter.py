import scipy.ndimage.filters as filters
import numpy as np
from skimage import io
import matplotlib.pyplot as plt

# Load and normalize image
img = img = io.imread('SP_Brain.png',as_gray=True)
img = img/np.max(img)

#####
# initializing the filter of size 5 by 5
size_filter = 5
# the filter is divided by size_filter^2 for normalization
mean_filter = np.ones((size_filter,size_filter))/np.power(size_filter,2)
# Gaussain filter
gaussian_filter= filters.gaussian_filter(img, sigma=2) #, order=0, output=None, mode='reflect', cval=0.0, truncate=4.0)
# performing convolution
mean_img_filtered = filters.convolve(img, mean_filter,mode='mirror',cval=0.0)
#####

# performing the median filter
img_filtered = filters.median_filter(img,size=5,mode='nearest',cval=0.0)

# Display results
fig = plt.figure(figsize=(8, 5))
plt.subplot(221)
plt.imshow(img, cmap=plt.cm.gray)
plt.title('Original image'), plt.axis('off')
plt.subplot(222)
plt.imshow(img_filtered, cmap=plt.cm.gray)
plt.title('Median-filtered image'), plt.axis('off')
plt.subplot(223)
plt.imshow(mean_img_filtered, cmap=plt.cm.gray)
plt.title('Mean-filtered image'), plt.axis('off')
plt.subplot(224)
plt.imshow(gaussian_filter, cmap=plt.cm.gray)
plt.title('Gaussian-filtered image'), plt.axis('off')