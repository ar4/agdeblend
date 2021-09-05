import numpy as np
import matplotlib.pyplot as plt

nx = 16
ny = 24
n_traces = nx * ny
n_times = 512
true_data = np.zeros([32, 40, n_times], np.float32)
in_data = np.zeros([32, 40, n_times], np.float32)
out_data = np.zeros([32, 40, n_times], np.float32)
true_data[:16, :24] = np.fromfile(
    "data_4_true_data_0.bin", np.float32
).reshape(nx, ny, n_times)
true_data[16:, 16:] = np.fromfile(
    "data_4_true_data_1.bin", np.float32
).reshape(nx, ny, n_times)
in_data[:16, :24] = np.fromfile(
    "data_4_blended_data_0.bin", np.float32
).reshape(nx, ny, n_times)
in_data[16:, 16:] = np.fromfile(
    "data_4_blended_data_1.bin", np.float32
).reshape(nx, ny, n_times)
out_data[:16, :24] = np.fromfile(
    "out/data_4_deblended_data_0.bin", np.float32
).reshape(nx, ny, n_times)
out_data[16:, 16:] = np.fromfile(
    "out/data_4_deblended_data_1.bin", np.float32
).reshape(nx, ny, n_times)

tidx = 256

_, ax = plt.subplots(1, 5, figsize=(15, 7), sharex=True, sharey=True)
vmin, vmax = np.percentile(true_data, [2, 98])
ax[0].imshow(
    true_data[:, :, tidx].T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
)
ax[0].set_title("True")
ax[0].set_ylabel("y")
ax[0].set_xlabel("x")
ax[1].imshow(
    in_data[:, :, tidx].T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
)
ax[1].set_title("In")
ax[2].imshow(
    out_data[:, :, tidx].T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
)
ax[2].set_title("Out")
ax[3].imshow(
    in_data[:, :, tidx].T - out_data[:, :, tidx].T,
    aspect="auto",
    cmap="gray",
    vmin=vmin,
    vmax=vmax,
)
ax[3].set_title("In-Out")
ax[4].imshow(
    out_data[:, :, tidx].T - true_data[:, :, tidx].T,
    aspect="auto",
    cmap="gray",
    vmin=vmin,
    vmax=vmax,
)
ax[4].set_title("Out-True")
plt.show()
