import sys
import numpy as np
import matplotlib.pyplot as plt


def plot(
    n_traces,
    n_times,
    data_num,
    in_name,
    out_name,
    plotname,
    out_n_times=None,
    true_name=None,
    fileidx=None,
    out_fileidx=None,
):
    if out_n_times is None:
        out_n_times = n_times
    if fileidx is None:
        fileidx_str = ""
    else:
        fileidx_str = "_" + fileidx
    in_data = np.fromfile(
        "data_{}_{}_data{}.bin".format(data_num, in_name, fileidx_str),
        np.float32,
    ).reshape(n_traces, n_times)
    if out_fileidx is not None:
        fileidx_str = "_" + out_fileidx
    out_data = np.fromfile(
        "out/data_{}_{}_data{}.bin".format(data_num, out_name, fileidx_str),
        np.float32,
    ).reshape(n_traces, out_n_times)
    if true_name is not None:
        true_data = np.fromfile(
            "data_{}_{}_data{}.bin".format(data_num, true_name, fileidx_str),
            np.float32,
        ).reshape(n_traces, n_times)
        n_plots = 3
    else:
        n_plots = 2

    _, ax = plt.subplots(
        1, n_plots, figsize=(7, 3.5), sharex=True, sharey=True
    )
    if out_n_times != n_times:
        in_data = np.pad(in_data, ((0, 0), (0, out_n_times - n_times)))
    vmin, vmax = np.percentile(in_data, [2, 98])
    ax[0].imshow(in_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
    ax[0].set_title("In")
    ax[0].set_ylabel("Time sample")
    ax[0].set_xlabel("Trace")
    ax[1].imshow(out_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax)
    ax[1].set_title("Out")
    if true_name is not None:
        ax[2].imshow(
            true_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
        )
        ax[2].set_title("True")
    plt.tight_layout()
    plt.savefig(plotname)


def plot_wavelet():
    wavelet = np.fromfile("data_2_wavelet.bin", np.float32)
    plt.figure(figsize=(7, 3.5))
    plt.plot(wavelet)
    plt.xlabel("Time sample")
    plt.ylabel("Amplitude")
    plt.tight_layout()
    plt.savefig("example_7_wavelet.png")


def plot_9_traces():
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

        _, ax = plt.subplots(1, 3, figsize=(7, 3.5), sharex=True, sharey=True)
        vmin, vmax = np.percentile(true_data, [2, 98])
        ax[0].imshow(
            in_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
        )
        ax[0].set_title("In")
        ax[0].set_ylabel("Time sample")
        ax[0].set_xlabel("Trace")
        ax[1].imshow(
            out_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
        )
        ax[1].set_title("Out")
        ax[2].imshow(
            true_data.T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
        )
        ax[2].set_title("True")
        plt.tight_layout()
        plt.savefig("example_9_{}.jpg".format(fidx))


def plot_9_time():
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

    _, ax = plt.subplots(1, 3, figsize=(7, 3.5), sharex=True, sharey=True)
    vmin, vmax = np.percentile(true_data, [2, 98])
    ax[0].imshow(
        in_data[:, :, tidx].T, aspect="auto", cmap="gray", vmin=vmin, vmax=vmax
    )
    ax[0].set_title("In")

    ax[0].set_ylabel("y")
    ax[0].set_xlabel("x")
    ax[1].imshow(
        out_data[:, :, tidx].T,
        aspect="auto",
        cmap="gray",
        vmin=vmin,
        vmax=vmax,
    )
    ax[1].set_title("Out")
    ax[2].imshow(
        true_data[:, :, tidx].T,
        aspect="auto",
        cmap="gray",
        vmin=vmin,
        vmax=vmax,
    )
    ax[2].set_title("True")
    plt.tight_layout()
    plt.savefig("example_9_time.jpg")


if sys.argv[1] == "1":
    plot(128, 512, "1", "true", "blended", "example_1.jpg")
elif sys.argv[1] == "2":
    plot(
        128,
        512,
        "1",
        "blended",
        "blended",
        "example_2.jpg",
        out_n_times=768,
        out_fileidx="2",
    )
elif sys.argv[1] == "3":
    plot(128, 512, "1", "true", "blended", "example_3.jpg")
elif sys.argv[1] in ["4", "5", "6"]:
    plot(
        128,
        512,
        "1",
        "blended",
        "deblended",
        "example_" + sys.argv[1] + ".jpg",
        true_name="true",
    )
elif sys.argv[1] == "7":
    plot(
        64,
        1024,
        "2",
        "blended",
        "deblended",
        "example_7.jpg",
        true_name="true",
    )
    plot_wavelet()
elif sys.argv[1] == "8":
    plot(
        128,
        512,
        "3",
        "blended",
        "deblended",
        "example_8_0.jpg",
        true_name="true",
        fileidx="0",
    )
    plot(
        128,
        768,
        "3",
        "blended",
        "deblended",
        "example_8_1.jpg",
        true_name="true",
        fileidx="1",
    )
elif sys.argv[1] == "9":
    plot_9_traces()
    plot_9_time()
else:
    raise RuntimeError("unknown example index")
