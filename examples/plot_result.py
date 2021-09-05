import numpy as np
import matplotlib.pyplot as plt

n_traces = 16 * 16  # 128
n_times = 512
# true_data = np.fromfile("true_data_1.bin", np.float32).reshape(n_traces, n_times)
# blended_data = np.fromfile("blended_data_1.bin", np.float32).reshape(n_traces, n_times)
# deblended_data = np.fromfile("deblended_data_1.bin", np.float32).reshape(n_traces, n_times)
true_data = np.fromfile("data_2_true_data_1.bin", np.float32).reshape(
    n_traces, n_times
)
blended_data = np.fromfile("data_2_blended_data_1.bin", np.float32).reshape(
    n_traces, n_times
)
# deblended_data = np.fromfile("deblended_data_1.bin", np.float32).reshape(n_traces, n_times)

_, ax = plt.subplots(1, 2, figsize=(15, 7), sharex=True, sharey=True)
vmin, vmax = np.percentile(true_data, [2, 98])
ax[0].imshow(true_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
ax[0].set_title("True")
ax[0].set_ylabel("Time sample")
ax[0].set_xlabel("Trace")
ax[1].imshow(blended_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
ax[1].set_title("In")
# ax[2].imshow(deblended_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
# ax[2].set_title("Out")
# ax[3].imshow(blended_data.T-deblended_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
# ax[3].set_title("In-Out")
# ax[4].imshow(deblended_data.T-true_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
# ax[4].set_title("Out-True")
plt.show()
