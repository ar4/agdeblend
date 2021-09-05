import numpy as np
import matplotlib.pyplot as plt

nx = 16
ny = 24
n_traces = nx * ny
n_times = 512
for fidx in range(2):
    true_data = np.fromfile(
        "data_4_true_data_{}.bin".format(fidx), np.float32
    ).reshape(n_traces, n_times)
    in_data = np.fromfile(
        "data_4_blended_data_{}.bin".format(fidx), np.float32
    ).reshape(n_traces, n_times)
    out_data = np.fromfile(
        "out/data_4_deblended_data_{}.bin".format(fidx), np.float32
    ).reshape(n_traces, n_times)

    _, ax = plt.subplots(1, 5, figsize=(15, 7), sharex=True, sharey=True)
    vmin, vmax = np.percentile(true_data, [2, 98])
    ax[0].imshow(true_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
    ax[0].set_title("True")
    ax[0].set_ylabel("Time sample")
    ax[0].set_xlabel("Trace")
    ax[1].imshow(in_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
    ax[1].set_title("In")
    ax[2].imshow(out_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
    ax[2].set_title("Out")
    ax[3].imshow(
        in_data.T - out_data.T,
        aspect="auto",
        cmap="gray",
        vmin=vmin,
        vmax=vmax,
    )
    ax[3].set_title("In-Out")
    ax[4].imshow(
        out_data.T - true_data.T,
        aspect="auto",
        cmap="gray",
        vmin=vmin,
        vmax=vmax,
    )
    ax[4].set_title("Out-True")
plt.show()
