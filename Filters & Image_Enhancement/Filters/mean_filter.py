import scipy.ndimage.filters as filters
import numpy as np
from skimage import io
import matplotlib.pyplot as plt

# Load and normalize image
img = img = io.imread('Brain.tif',as_gray=True)
img = img/np.max(img)

# initializing the filter of size 5 by 5
size_filter = 5
# the filter is divided by size_filter^2 for normalization
mean_filter = np.ones((size_filter,size_filter))/np.power(size_filter,2)

# Gaussain filter
gaussian_filter= filters.gaussian_filter(img, sigma=2) #, order=0, output=None, mode='reflect', cval=0.0, truncate=4.0)

# performing convolution
img_filtered = filters.convolve(img, mean_filter,mode='constant',cval=0.0)

# Display results
fig = plt.figure(figsize=(8, 5))
plt.subplot(131)
plt.imshow(img, cmap=plt.cm.gray)
plt.title('Original image'), plt.axis('off')
plt.subplot(132)
plt.imshow(img_filtered, cmap=plt.cm.gray)
plt.title('Mean-filtered image'), plt.axis('off')
plt.subplot(133)
plt.imshow(gaussian_filter, cmap=plt.cm.gray)
plt.title('Gaussian-filtered image'), plt.axis('off')