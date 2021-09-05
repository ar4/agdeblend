import numpy as np
import segyio
import matplotlib.pyplot as plt
from agdeblend import blend, AGDBlendSum

n_shots = 128
n_channels = 120
n_times = 1500

with segyio.open(
    "mobil_avo_viking_graben_line_12.segy", ignore_geometry=True
) as s:
    d_true = s.trace.raw[: n_shots * n_channels].reshape(
        n_shots, n_channels, n_times
    )

delay = 750
jitter = 50
shottimes = (
    np.arange(n_shots) * delay + (np.random.rand(n_shots) - 0.5) * 2 * jitter
).astype(int)
shottimes = np.tile(shottimes.reshape(-1, 1), [1, n_channels])
channels = np.tile(np.arange(n_channels).reshape(1, -1), [n_shots, 1])
trace_types = np.zeros([n_shots, n_channels])
blend_mode = AGDBlendSum
d_blended = blend(
    [shottimes], [channels], [trace_types], [d_true], blend_mode
)[0]

# Blended data in the shot and channel domains
_, ax = plt.subplots(1, 2, figsize=(7, 3.5), sharey=True)
vmin, vmax = np.percentile(d_blended[10], [5, 95])
ax[0].imshow(
    d_blended[10].T,
    aspect="auto",
    cmap="gray",
    vmin=vmin,
    vmax=vmax,
    interpolation="bilinear",
)
ax[0].set_xlabel("Channel index")
ax[0].set_ylabel("Time sample")
ax[0].set_title("Shot gather")
vmin, vmax = np.percentile(d_blended[:, -1], [5, 95])
ax[1].imshow(
    d_blended[:, -1].T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
)
ax[1].set_xlabel("Shot index")
ax[1].set_title("Channel gather")
plt.tight_layout()
plt.savefig("blended.jpg")

# Signal is sparse in the Fourier domain
D_TRUE = np.abs(np.fft.rfft2(d_true[:64, -1]))
D_BLENDED = np.abs(np.fft.rfft2(d_blended[:64, -1]))
_, ax = plt.subplots(2, 2, figsize=(7, 7))
vmin, vmax = np.percentile(d_true[:64, -1], [5, 95])
ax[0, 0].imshow(
    d_true[:64, -1].T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
)
ax[0, 0].set_xlabel("Shot index")
ax[0, 0].set_ylabel("Time index")
ax[0, 0].set_title("Signal TX")
ax[0, 1].imshow(
    d_blended[:64, -1].T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
)
ax[0, 1].set_title("Blended TX")
vmin, vmax = np.percentile(D_BLENDED, [5, 95])
ax[1, 0].imshow(D_TRUE.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
ax[1, 0].set_xlabel("Wavenumber index")
ax[1, 0].set_ylabel("Frequency index")
ax[1, 0].set_title("Signal FK")
ax[1, 1].imshow(D_BLENDED.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
ax[1, 1].set_title("Blended FK")
plt.tight_layout()
plt.savefig("transformed.jpg")

# Deblending is underdetermined
d1 = np.zeros([120, 1500 + 1500 // 2], np.float32)
d2 = np.zeros([120, 1500 + 1500 // 2], np.float32)
d1[:, :1500] = d_true[0]
d2[:, 1500 // 2 :] = d_true[1]
d_blended = d1 + d2
d_pd1 = np.zeros_like(d_blended)
d_pd2 = np.zeros_like(d_blended)
d_pd1[:, : 1500 - 375] = d_blended[:, : 1500 - 375]
d_pd2[:, 375:1500] = d_blended[:, 1500 - 375 :]

_, ax = plt.subplots(1, 3, figsize=(7, 3.5), sharex=True, sharey=True)
vmin, vmax = np.percentile(d_true[0], [10, 90])

ax[0].set_title("Shot 1")
ax[1].set_title("Shot 2")
ax[2].set_title("Recorded")
ax[0].set_yticks([])
ax[0].set_xticks([])
ax[0].set_ylabel("Time")
ax[0].set_xlabel("Channel")
ax[0].imshow(
    d_pd1.T,
    aspect="auto",
    cmap="gray",
    vmin=vmin,
    vmax=vmax,
    interpolation="bilinear",
)
ax[1].imshow(
    d_pd2.T,
    aspect="auto",
    cmap="gray",
    vmin=vmin,
    vmax=vmax,
    interpolation="bilinear",
)
ax[2].imshow(
    d_blended.T,
    aspect="auto",
    cmap="gray",
    vmin=vmin,
    vmax=vmax,
    interpolation="bilinear",
)
plt.tight_layout()
plt.savefig("underdetermined.jpg")

# Example 9 patches


def add_patch(patch=None):
    volume = np.zeros([32, 40], int)
    file_0_range = (slice(None, 16), slice(None, 24))
    file_1_range = (slice(16, 32), slice(16, 40))
    volume[file_0_range] = 1
    volume[file_1_range] = 1
    if patch is not None:
        volume[patch] += 2
    return volume.T


patch_ranges_x = [slice(0, 16), slice(8, 24), slice(16, 32)]
patch_ranges_y = [slice(0, 16), slice(8, 32), slice(24, 40)]
patch_ranges = []
for ix in range(3):
    for iy in range(3):
        if (ix == 2 and iy == 0) or (ix == 0 and iy == 2):
            continue
        patch_ranges.append((patch_ranges_x[ix], patch_ranges_y[iy]))

_, ax = plt.subplots(2, 4, figsize=(7, 7), sharex=True, sharey=True)
vmin = 0
vmax = 3
ax[0, 0].imshow(add_patch(), vmin=vmin, vmax=vmax)
ax[0, 0].set_title("Volume")
for i, patch_range in enumerate(patch_ranges):
    j = i + 1
    row = j // 4
    col = j % 4
    ax[row, col].imshow(add_patch(patch_range), vmin=vmin, vmax=vmax)
    ax[row, col].set_xticks([])
    ax[row, col].set_yticks([])
    ax[row, col].set_title("Patch {}".format(j))

plt.tight_layout()
plt.savefig("example_9_patches.png")
