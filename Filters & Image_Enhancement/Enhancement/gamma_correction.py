import numpy as np
from skimage import io
import plot_tools as plt_hist
import matplotlib.pyplot as plt

# Load and normalize image
img = img = io.imread('Chest_RX_bright.png',as_gray=True)
img = img/np.max(img)

# gamma is initialized
gamma = 2

# gamma-correction exponent is computed
img_log = np.log(img)*gamma

# gamma-correction is performed
img_gamma_corrected = np.exp(img_log)

# Display results
fig = plt.figure(figsize=(8, 5))
axes = np.zeros((2, 2), dtype=np.object)
axes[0, 0] = fig.add_subplot(2, 2, 1)
for i in range(1, 2):
    axes[0, i] = fig.add_subplot(2, 2, 1+i, sharex=axes[0,0], sharey=axes[0,0])
for i in range(0, 2):
    axes[1, i] = fig.add_subplot(2, 2, 3+i)
    
ax_img, ax_hist = plt_hist.plot_img_and_hist(img, axes[:, 0])
ax_img.set_title('Low contrast bright image')

y_min, y_max = ax_hist.get_ylim()
ax_hist.set_ylabel('Number of pixels')
ax_hist.set_yticks(np.linspace(0, y_max, 5))

ax_img, ax_hist = plt_hist.plot_img_and_hist(img_gamma_corrected, axes[:, 1])
ax_img.set_title('Gamma correction')

# prevent overlap of y-axis labels
fig.tight_layout()
plt.show()