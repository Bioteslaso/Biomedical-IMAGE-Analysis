import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from skimage import exposure
from skimage import io
import plot_tools as plt_hist

matplotlib.rcParams['font.size'] = 8

# Load the image
img = io.imread('Brain.tif',as_gray=True)

# Since the result of the equalization is a normalized image, we also normalize
# the original image.
img = img/np.max(img)

# Equalization
img_eq = exposure.equalize_hist(img)

# Display results
fig = plt.figure(figsize=(8, 5))
axes = np.zeros((2, 2), dtype=np.object)
axes[0, 0] = fig.add_subplot(2, 2, 1)
for i in range(1, 2):
    axes[0, i] = fig.add_subplot(2, 2, 1+i, sharex=axes[0,0], sharey=axes[0,0])
for i in range(0, 2):
    axes[1, i] = fig.add_subplot(2, 2, 3+i)

ax_img, ax_hist, ax_cdf = plt_hist.plot_img_and_hist(img, axes[:, 0], printCDF=True)
ax_img.set_title('Low contrast image')

y_min, y_max = ax_hist.get_ylim()
ax_hist.set_ylabel('Number of pixels')
ax_hist.set_yticks(np.linspace(0, y_max, 5))

ax_img, ax_hist, ax_cdf = plt_hist.plot_img_and_hist(img_eq, axes[:, 1], printCDF=True)
ax_img.set_title('Histogram equalization')

ax_cdf.set_ylabel('Fraction of total intensity')
ax_cdf.set_yticks(np.linspace(0, 1, 5))

# prevent overlap of y-axis labels
fig.tight_layout()
plt.show()