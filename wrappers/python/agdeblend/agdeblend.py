"""AGDeblend seismic data blending and deblending."""
import ctypes
import numpy as np

AGDLive = 0
AGDBad = 1
AGDMissing = 2
AGDBlendSum = 0
AGDBlendMean = 1
AGDBlendOverwrite = 2


def deblend(
    volumes,
    window_shapes,
    coords,
    shottimes,
    channels,
    trace_types,
    initial_factor,
    n_its,
    print_freq,
    data,
    wavelet_idxs=None,
    wavelets=None,
    comm=None,
):
    """Deblends seismic data.

    See https://ausargeo.com/agdeblend for documentation.

    Args:
      volumes: [n_patches] int array specifying the volume index that each patch
          belongs to
      window_shapes: [n_volume][<number of space dimensions in volume>] in array
          specifying the window shapes to use
      coords: [n_patches][<number of space dimensions in patch>] int array
          specifying the coordinates of each patch within its volume
      shottimes: [n_patches][<space shape of patch>] long int array specifying
          the time sample index of each trace's shot in the continuous record
      channels: [n_patches][<space shape of patch>] int array specifying the
          recording channel number of each trace
      trace_types: [n_patches][<space shape of patch>] int array specifying
          if the trace is live (AGDLive), bad (AGDBad) or missing (AGDMissing)
      initial_factor: Float specifying the initial threshold factor in (0, 1]
      n_its: Int specifying the number of iterations
      print_freq: Int specifying the iteration interval between
          printing the norm of the residual, with values less than 1
          corresponding to never
      data: [n_patches][<shape of patch>] blended seismic data (may be
          modified during deblending)
      wavelet_idxs: [n_patches][<space shape of patch>] optional int array
          specifying the index in the wavelets array (see below) to use
          as the source wavelet for each trace (default None)
      wavelets: [number of wavelets][number of time samples in each wavelet]
          optional array containing time samples of source wavelets
          (default None)
      comm: Optional MPI Communicator (default None, which will run a
          serial implementation)

      Returns:
        The deblended data in an array of the same size as the input
    """

    fdtype, cfdtype, fdtype_ptrptr = _get_fdtype(data[0])
    agd_lib, ffi, comm = _get_agd_lib(cfdtype, comm)

    # Convert to contiguous NumPy arrays of the right type
    n_patches = len(data)

    data_shapes = _convert_to_numpy([patch_data.shape for patch_data in data])
    data = _convert_to_numpy(data, dtype=cfdtype, extra_requirements=["W"])
    volumes = np.require(volumes, dtype=ctypes.c_int, requirements=["C", "A"])
    coords = _convert_to_numpy(coords)
    window_shapes = _convert_to_numpy(window_shapes)
    shottimes = _convert_to_numpy(shottimes, dtype=ctypes.c_long)
    channels = _convert_to_numpy(channels)
    trace_types = _convert_to_numpy(trace_types)
    if wavelet_idxs is not None or wavelets is not None:
        if wavelet_idxs is None or wavelets is None:
            raise RuntimeError(
                "If wavelet_idxs or wavelets is not None, they both must not be"
            )
        wavelet_lengths = np.array(
            [len(wavelet) for wavelet in wavelets]
        ).astype(ctypes.c_int)
        wavelet_lengths_ptr = wavelet_lengths.ctypes.data
        wavelet_idxs = _convert_to_numpy(wavelet_idxs)
        wavelet_idxs_ptr = wavelet_idxs["ppointers"]
        wavelets = _convert_to_numpy(wavelets, dtype=cfdtype)
        wavelets_ptr = wavelets["ppointers"]
    else:
        wavelet_lengths_ptr = ffi.NULL
        wavelet_idxs_ptr = ffi.NULL
        wavelets_ptr = ffi.NULL

    # Check that the number of patches is the same for all variables
    if len(volumes) != n_patches:
        raise RuntimeError("volumes must have the same length as data")
    _check_n_patches(n_patches, coords, "coords")
    _check_n_patches(n_patches, shottimes, "shottimes")
    _check_n_patches(n_patches, channels, "channels")
    _check_n_patches(n_patches, trace_types, "trace_types")
    if wavelet_idxs is not None:
        _check_n_patches(n_patches, wavelet_idxs, "wavelet_idxs")

    # Check that the number of traces in each patch is the same for
    # all variables
    n_traces = [np.prod(patch_data.shape[:-1]) for patch_data in data["np"]]
    _check_n_traces(n_traces, shottimes, "shottimes")
    _check_n_traces(n_traces, channels, "channels")
    _check_n_traces(n_traces, trace_types, "trace_types")
    if wavelet_idxs is not None:
        _check_n_traces(n_traces, wavelet_idxs, "wavelet_idxs")

    # Check wavelet_idxs match wavelets
    if wavelet_idxs is not None:
        _check_wavelets(wavelet_idxs, wavelets, trace_types)

    # Extract the number of dimensions in each volume
    n_volumes = np.max(volumes) + 1
    n_dims = np.zeros(n_volumes, ctypes.c_int)
    for volume_idx, patch_data in zip(volumes, data["np"]):
      n_dims[volume_idx] = patch_data.ndim

    if agd_lib.agd_deblend(
        n_patches,
        ffi.cast("int *", volumes.ctypes.data),
        ffi.cast("int *", n_dims.ctypes.data),
        ffi.cast("int **", window_shapes["ppointers"]),
        ffi.cast("int **", coords["ppointers"]),
        ffi.cast("int **", data_shapes["ppointers"]),
        ffi.cast("long int **", shottimes["ppointers"]),
        ffi.cast("int **", channels["ppointers"]),
        ffi.cast("enum AGDTraceType **", trace_types["ppointers"]),
        ffi.cast("int *", wavelet_lengths_ptr),
        ffi.cast("int **", wavelet_idxs_ptr),
        ffi.cast(fdtype_ptrptr, wavelets_ptr),
        initial_factor,
        n_its,
        print_freq,
        *comm,
        ffi.cast(fdtype_ptrptr, data["ppointers"])
    ):
        raise RuntimeError("Error in call to agd_deblend")

    return [patch_data.astype(fdtype) for patch_data in data["np"]]


def blend(
    shottimes,
    channels,
    trace_types,
    data,
    blend_mode,
    taper_length=0,
    n_times_out=None,
    shottimes_out=None,
    channels_out=None,
    trace_types_out=None,
    comm=None,
    data_out=None,
):
    """Blends seismic data or modifies already blended data.

    See https://ausargeo.com/agdeblend for documentation.

    Args:
      shottimes: [n_patches][<space shape of patch>] long int array specifying
          the time sample index of each trace's shot in the continuous record
      channels: [n_patches][<space shape of patch>] int array specifying the
          recording channel number of each trace
      trace_types: [n_patches][<space shape of patch>] int array specifying
          if the trace is live (AGDLive), bad (AGDBad) or missing (AGDMissing)
      blend_mode: int specifying whether to sum overlapping samples
          (AGDBlendSum), overwrite them (AGDBlendOverwrite), or take the mean
          (AGDBlendMean)
      taper_length: Optional int specifying the taper length to use in time
          samples at the beginning and end of each trace if using the
          AGDBlendMean blend mode (default 0)
      data: [n_patches][<shape of patch>] seismic data (blended or unblended)
      n_times_out: [n_patches] optional int array specifying the number of
          time samples to extract from the continuous record after blending
          for each output patch (defaults to the same as the input)
      shottimes_out: [n_patches_out][<space shape of patch>] optional long int
          array specifying the time samples indexes of the continuous record
          to extract each output trace from after blending (defaults to the
          same as the input)
      channels_out: [n_patches_out][<space shape of patch>] optional int array
          specifying the channel to extract each output trace from after
          blending (defaults to the same as the input)
      trace_types_out: [n_patches_out][<space shape of patch>] optional int
          array specifying the trace type of each output trace, with
          AGDLive and AGDMissing being the only allowed options (defaults
          to the same as the input)
      comm: Optional MPI Communicator (default None, which will run a
          serial implementation)
      data_out: [n_patches_out][<shape of output patch>] optional output
          array, which may be the input data array if the output will
          be the same shape (default None, which will allocate a new array)

      Returns:
        The data after blending to a continuous record and extracting the output
        traces as requested, in an array of size
        [n_patches_out][<shape of output patch>]
    """

    fdtype, cfdtype, fdtype_ptrptr = _get_fdtype(data[0])
    agd_lib, ffi, comm = _get_agd_lib(cfdtype, comm)

    if blend_mode not in [AGDBlendSum, AGDBlendMean, AGDBlendOverwrite]:
        raise RuntimeError("Invalid blend_mode")

    # Convert to contiguous NumPy arrays of the right type
    data = _convert_to_numpy(data, dtype=cfdtype)
    shottimes = _convert_to_numpy(shottimes, dtype=ctypes.c_long)
    channels = _convert_to_numpy(channels)
    trace_types = _convert_to_numpy(trace_types)
    shottimes_out = _convert_to_numpy(
        shottimes_out, dtype=ctypes.c_long, default=shottimes
    )
    channels_out = _convert_to_numpy(channels_out, default=channels)
    trace_types_out = _convert_to_numpy(trace_types_out, default=trace_types)

    # Set n_patches, n_traces, n_times, + _out versions, and allocate out_data
    n_patches = len(data["np"])
    n_traces = np.zeros(n_patches, ctypes.c_int)
    n_times = np.zeros(n_patches, ctypes.c_int)
    for i, patch_data in enumerate(data["np"]):
        n_traces[i] = np.prod(patch_data.shape[:-1])
        n_times[i] = patch_data.shape[-1]

    if data_out is not None:
        if len(data_out) < 1:
            raise RuntimeError("data_out must have at least one patch")
        if not isinstance(data_out[0], np.ndarray):
            raise RuntimeError("data_out patches must be NumPy arrays")
        if data_out[0].dtype != cfdtype:
            raise RuntimeError("data_out does not have the right data type")
        n_patches_out = len(data_out)
        n_traces_out = np.array(
            [np.prod(patch_data.shape[:-1]) for patch_data in data_out]
        ).astype(ctypes.c_int)
        if n_times_out is None:
            n_times_out = [patch_data.shape[-1] for patch_data in data_out]
        else:
            if not all(
                n_times_out
                == [patch_data.shape[-1] for patch_data in data_out]
            ):
                raise RuntimeError("data_out does not match n_times_out")
        n_times_out = np.array(n_times_out).astype(ctypes.c_int)
        data_out = {"np": data_out}
    else:  # data_out not provided, so use shottimes_out and n_times_out/n_times
        n_patches_out = len(shottimes_out["np"])
        n_traces_out = np.array(
            [
                np.prod(patch_shottimes.shape)
                for patch_shottimes in shottimes_out["np"]
            ]
        ).astype(ctypes.c_int)
        if n_times_out is None:
            n_times_out = n_times
        else:
            n_times_out = np.array(n_times_out).astype(ctypes.c_int)
        data_out = {
            "np": [
                np.zeros(
                    patch_shottimes.shape + (patch_n_times,), dtype=cfdtype
                )
                for patch_shottimes, patch_n_times in zip(
                    shottimes_out["np"], n_times_out
                )
            ]
        }

    if len(n_times_out) != n_patches_out:
        raise RuntimeError("n_times_out is not of the expected length")

    data_out["pointers"] = np.array(
        [patch_var.ctypes.data for patch_var in data_out["np"]]
    )
    data_out["ppointers"] = data_out["pointers"].ctypes.data

    # Check that the number of patches is the same for all variables
    _check_n_patches(n_patches, shottimes, "shottimes")
    _check_n_patches(n_patches, channels, "channels")
    _check_n_patches(n_patches, trace_types, "trace_types")
    _check_n_patches(n_patches_out, shottimes_out, "shottimes_out")
    _check_n_patches(n_patches_out, channels_out, "channels_out")
    _check_n_patches(n_patches_out, trace_types_out, "trace_types_out")

    # Check that the number of traces in each patch is the same for
    # all variables
    _check_n_traces(n_traces, shottimes, "shottimes")
    _check_n_traces(n_traces, channels, "channels")
    _check_n_traces(n_traces, trace_types, "trace_types")
    _check_n_traces(n_traces_out, shottimes_out, "shottimes_out")
    _check_n_traces(n_traces_out, channels_out, "channels_out")
    _check_n_traces(n_traces_out, trace_types_out, "trace_types_out")

    if agd_lib.agd_blend(
        n_patches,
        ffi.cast("int *", n_traces.ctypes.data),
        ffi.cast("int *", n_times.ctypes.data),
        ffi.cast("long int **", shottimes["ppointers"]),
        ffi.cast("int **", channels["ppointers"]),
        ffi.cast("enum AGDTraceType **", trace_types["ppointers"]),
        ffi.cast(fdtype_ptrptr, data["ppointers"]),
        ffi.cast("enum AGDBlendMode", blend_mode),
        taper_length,
        n_patches_out,
        ffi.cast("int *", n_traces_out.ctypes.data),
        ffi.cast("int *", n_times_out.ctypes.data),
        ffi.cast("long int **", shottimes_out["ppointers"]),
        ffi.cast("int **", channels_out["ppointers"]),
        ffi.cast("enum AGDTraceType **", trace_types_out["ppointers"]),
        *comm,
        ffi.cast(fdtype_ptrptr, data_out["ppointers"])
    ):
        raise RuntimeError("Error in call to agd_blend")

    return [patch_data.astype(fdtype) for patch_data in data_out["np"]]


def _get_fdtype(np_var):
    fdtype = np_var.dtype
    if fdtype not in (np.float32, np.float64):
        raise RuntimeError("Data must be float32 or float64")
    if fdtype == np.float32:
        cfdtype = ctypes.c_float
        fdtype_ptrptr = "float **"
    else:
        cfdtype = ctypes.c_double
        fdtype_ptrptr = "double **"
    return fdtype, cfdtype, fdtype_ptrptr


def _get_agd_lib(cfdtype, comm):
    if cfdtype == ctypes.c_float and comm is None:
        try:
            import agdeblend._agdeblend_ffi

            agd_lib = agdeblend._agdeblend_ffi.lib
            ffi = agdeblend._agdeblend_ffi.ffi
            comm = []
        except:
            raise RuntimeError(
                "The single precision, non-MPI AGDeblend library was not found"
            )
    elif cfdtype == ctypes.c_double and comm is None:
        try:
            import agdeblend._agdeblend_double_ffi

            agd_lib = agdeblend._agdeblend_double_ffi.lib
            ffi = agdeblend._agdeblend_double_ffi.ffi
            comm = []
        except:
            raise RuntimeError(
                "The double precision, non-MPI AGDeblend library was not found"
            )
    elif cfdtype == ctypes.c_float and comm is not None:
        try:
            import agdeblend._agdeblend_mpi_ffi
            from mpi4py import MPI

            agd_lib = agdeblend._agdeblend_mpi_ffi.lib
            ffi = agdeblend._agdeblend_mpi_ffi.ffi
            comm = [ffi.cast("MPI_Comm *", MPI._addressof(comm))[0]]
        except:
            raise RuntimeError(
                "The single precision, MPI AGDeblend library was not found"
            )
    elif cfdtype == ctypes.c_double and comm is not None:
        try:
            import agdeblend._agdeblend_mpi_double_ffi
            from mpi4py import MPI

            agd_lib = agdeblend._agdeblend_mpi_double_ffi.lib
            ffi = agdeblend._agdeblend_mpi_double_ffi.ffi
            comm = [ffi.cast("MPI_Comm *", MPI._addressof(comm))[0]]
        except:
            raise RuntimeError(
                "The double precision, MPI AGDeblend library was not found"
            )

    return agd_lib, ffi, comm


def _convert_to_numpy(var, dtype=None, extra_requirements=None, default=None):
    if var is None:
        return default
    np_var = []
    if dtype is None:
        dtype = ctypes.c_int
    requirements = ["C_CONTIGUOUS", "ALIGNED"]
    if extra_requirements is not None:
        requirements += extra_requirements
    for patch_var in var:
        np_var.append(
            np.require(patch_var, dtype=dtype, requirements=requirements)
        )
    pointers = np.array([patch_var.ctypes.data for patch_var in np_var])
    return {
        "np": np_var,
        "pointers": pointers,
        "ppointers": pointers.ctypes.data,
    }


def _check_n_patches(n_patches, var, name):
    if len(var["np"]) != n_patches:
        raise RuntimeError(name + " must have the same length as data")


def _get_n_traces(var):
    return [np.prod(patch_var.shape) for patch_var in var["np"]]


def _check_n_traces(n_traces, var, name):
    if _get_n_traces(var) != n_traces:
        raise RuntimeError(name + " has an incorrect number of traces")


def _check_wavelets(wavelet_idxs, wavelets, trace_types):
    min_wavelet_idx = 0
    max_wavelet_idx = 0
    for patch_wavelet_idxs, patch_trace_types in zip(
        wavelet_idxs["np"], trace_types["np"]
    ):
        for trace_wavelet_idx, trace_type in zip(
            patch_wavelet_idxs.ravel(), patch_trace_types.ravel()
        ):
            if trace_type == AGDLive:
                min_wavelet_idx = min(trace_wavelet_idx, min_wavelet_idx)
                max_wavelet_idx = min(trace_wavelet_idx, max_wavelet_idx)
    if min_wavelet_idx < 0 or max_wavelet_idx >= len(wavelets["np"]):
        raise RuntimeError("wavelet_idx must be in [0, len(wavelets))")
