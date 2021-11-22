Usage
=====

AGDeblend is designed to be easy to use for simple cases, but to provide flexibility for more advanced cases, which can be a bit more complicated. The :doc:`examples` demonstrate usage of AGDeblend for a range of situations in C, Python, Julia, and Fortran. If you haven't done so already, you should also read the :doc:`introduction` to learn the terminology and concepts that AGDeblend uses.

Blend
-----

The `agd_blend` function takes input traces, blends them into continuous records, and then pseudo-deblends them back into traces.

It can be used to create synthetic data by using unblended traces as the input and the :c:enumerator:`AGDBlendSum` blend mode, which will sum overlapping traces, as in :doc:`example_1`.

It can also be used to adjust already blended data. For example, if you wish to extend the length of traces, you could use the blended data and the input, with the :c:enumerator:`AGDBlendMean` or :c:enumerator:`AGDBlendOverwrite` blend mode, and set output trace lengths (`n_times_out`) to be the desired length. This should recreate the continuous record and extract longer traces from it for the output. This is what is performed in :doc:`example_2`.

Another example adjustment of existing blended data is to extract traces at additional times, such as if there was interference from a nearby survey and you wish to extract traces at the shot times of that survey so that they can be included in deblending. In that case you would add the new shot times to `shottimes_out`.

Finally, this function is also useful to ensure that samples in the pseudo-deblended data that should be duplicates of each other, as they come from the same sample in the continuous record, really are equal. It is a good idea to do this before deblending if you have applied any processing to the blended data.

.. c:function:: int agd_blend(int n_patches, int const *n_traces, int const *n_times, long int const *const *shottimes, int const *const *channels, enum AGDTraceType const *const *trace_types, AGD_TYPE *const *data, enum AGDBlendMode blend_mode, int taper_length, int n_patches_out, int const *n_traces_out, int const *n_times_out, long int const *const *shottimes_out, int const *const *channels_out, enum AGDTraceType const *const *trace_types_out, MPI_Comm comm, AGD_TYPE *const *data_out)

  :param n_patches: The number of input patches
  :param n_traces: The number of traces in each input patch
  :param n_times: The number of samples in the time dimension of each input patch
  :param shottimes: The shot time (in units of samples) for each input trace
  :param channels: The channel label for each input trace
  :param trace_types: The trace type (:c:enumerator:`AGDLive`, :c:enumerator:`AGDBad`, or :c:enumerator:`AGDMissing`) for each input trace
  :param data: The input data
  :param blend_mode: The blend mode: :c:enumerator:`AGDBlendSum`, :c:enumerator:`AGDBlendMean`, or :c:enumerator:`AGDBlendOverwrite`
  :param taper_length: Only used with blend mode :c:enumerator:`AGDBlendMean`, when it determines how long the taper is (in units of time samples from the nearest end of the trace) of the weight used for calculating the weighted mean over overlapping samples (may be 0)
  :param n_patches_out: The number of output patches
  :param n_traces_out: The number of traces in each output patch
  :param n_times_out: The number of samples in the time dimension of each output patch
  :param shottimes_out: The shot time (in units of samples) for each output trace
  :param channels_out: The channel label for each output trace
  :param trace_types_out: The trace type (:c:enumerator:`AGDLive` or :c:enumerator:`AGDMissing`, but not :c:enumerator:`AGDBad`) for each output trace
  :param comm: The MPI Communicator to use (this parameter is only present if compiled with MPI support)
  :param data_out: The output data
  :returns: Zero on success, non-zero on failure

Arguments that are supposed to contain a value for each patch, such as `n_traces`, should be arrays of length `n_patches`. For the parameters that provide a value for each trace, a pointer to an array of pointers of length `n_patches` should be used. Each of the pointers in this array should point to arrays with a number of elements equal to the number of traces in the corresponding patch. Similarly, the input and output data arguments should be pointers to an array of `n_patches` pointers, each pointing to memory containing the data samples for that patch.

The same pointer can be used for the input and output arrays if they are the same, for example if the input and output shot times are the same then you can just pass the same pointer for both. If the memory pointed to by the input data is sufficient to also store the output data, then the same pointer may be used for both of those as well (and the input data will be overwritten).

`AGD_TYPE` will be replaced by `double` if the `AGD_DOUBLE` macro is defined during compilation, and `float` otherwise.

The origin of the times provided by `shottimes` is arbitrary - only the relative times are used.

The channel labels should be a unique number for each channel, but do not need to be ordered, sequential, or positive. 

.. c:enum:: AGDTraceType

  .. c:enumerator:: AGDLive
    
    Regular live trace
      
  .. c:enumerator:: AGDBad

    A bad trace that will not be used, and any samples that it overlaps with will be muted

  .. c:enumerator:: AGDMissing

    A missing trace that will be ignored when blending into a continuous record

.. c:enum:: AGDBlendMode

  .. c:enumerator:: AGDBlendSum

    Overlapping samples will be summed, as in real blending

  .. c:enumerator:: AGDBlendMean

    The weighted mean over overlapping samples will be used, with the weight determined by the distance from the nearest end of the trace and the taper length

  .. c:enumerator:: AGDBlendOverwrite

    Overlapping samples will overwrite each other, so the final value will be the last one written

Deblend
-------

The interface of the deblending function, `agd_deblend`, is similar to that of the blending one. 

.. c:function:: int agd_deblend(int n_patches, int const *volumes, int const *n_dims, int const *const *window_shapes, int const *const *coords, int const *const *shapes, long int const *const *shottimes, int const *const *channels, enum AGDTraceType const *const *trace_types, int const *wavelet_lengths, int const *const *wavelet_idxs, AGD_TYPE const *const *wavelets, AGD_TYPE const initial_factor, int n_its, int print_freq, MPI_Comm comm, AGD_TYPE *const *data)

  :param n_patches: The number of input patches
  :param volumes: The index of the volume that each patch belongs to (sequentially increasing from 0)
  :param n_dims: The number of dimensions in each volume
  :param window_shapes: The shapes of windows used for each volume
  :param coords: The coordinates of each patch within its volume
  :param shapes: The shapes (number of elements in each dimension) of each patch
  :param shottimes: The shot time (in units of samples) for each trace
  :param channels: The channel label for each trace
  :param trace_types: The trace type (:c:enumerator:`AGDLive`, :c:enumerator:`AGDBad`, or :c:enumerator:`AGDMissing`) for each trace
  :param wavelet_lengths: The number of time samples in each wavelet (NULL if no wavelets)
  :param wavelet_idxs: The index within the wavelet array of the wavelet to use for each trace (NULL if no wavelets)
  :param wavelets: An array of source wavelets to convolve with the traces prior to blending (NULL if no wavelets)
  :param initial_factor: An initial factor, in the range :math:`(0, 1]`, to apply to the threshold value
  :param n_its: The number of iterations
  :param print_freq: The number of iterations between printing the norm of the residual (0 or below means never)
  :param comm: The MPI Communicator to use (this parameter is only present if compiled with MPI support)
  :param data: The pseudo-deblended input data, and also where the deblended output will be written
  :returns: Zero on success, non-zero on failure

Unlike `coords` and `shottimes`, which have arbitrary origins, the volume indexes provided in `volumes` must correspond to sequentially numbered volumes beginning at 0. This is so that the values in `n_dims` and `window_shapes` can be associated with the correct volume. For example, if you only have one volume, the value in `volumes` for every patch should be 0, `n_dims` should contain one value, and `window_shapes` should contain a pointer to one array, of length equal to the value in `n_dims`.

When MPI is used, the volume indexes must be consistent across processes: for example volume 0 must correspond to the same volume on every process (even if a process does not contain any patches in that volume). This also means that the `window_shapes` and `n_dims` arrays should be the same on all processes. Wavelet indexes, however, do not need to be consistent. A trace in the overlap region between patches assigned to two different processes may have a different wavelet index on each, as long as the actual wavelet that those indexes correspond to on each process is the same.

More iterations will generally produce a better result. The appropriate number depends on the features of each dataset, such as how well separated the signal and blending noise are in the Fourier domain, and the blending factor. For a typical deep marine tower streamer with low noise and a blending factor of three, five hundred iterations may be sufficient, while in more difficult situations one thousand or more might be required.

Wrappers
--------

The interfaces provided by the wrappers for different programming languages are quite similar to the C interface, but there are some differences that are documented below.

The source code for the wrappers can be found in the `wrappers` directory. All of the :doc:`examples` include Python, Julia, and Fortran implementations, which should hopefully help to clarify usage.

Python
^^^^^^
.. py:function:: blend(shottimes, channels, trace_types, data, blend_mode, taper_length=0, n_times_out=None, shottimes_out=None, channels_out=None, trace_types_out=None, comm=None, data_out=None)

  :returns: The blended data

The parameters have the same meanings as their equivalents in the C interface. A value of `None` for the output parameters signifies that you do not wish to change it from the corresponding input argument. The available blend modes and trace types can be imported from the module.

.. py:function:: deblend(volumes, window_shapes, coords, shottimes, channels, trace_types, initial_factor, n_its, print_freq, data, wavelet_idxs=None, wavelets=None, comm=None)

  :returns: The deblended data

If possible, the input data memory will be reused to store the output, so the contents of the input data may be modified.

Julia
^^^^^

The Julia interface is the same as that of the Python wrapper.

As Julia is a column-major language, some arrays (such as the window shapes) will be specified the other way around compared to C and Python. The wrapper takes care of flipping these back to the way expected by the C code. One part of the interface that may feel somewhat unnatural in Julia, however, is that the volume indexes should still be zero-indexed, as in the C code.

Fortran
^^^^^^^

To fit more naturally with typical Fortran style, the Fortran interface is somewhat different from the other languages.

.. f:subroutine:: blend(patches, blend_mode, comm, ierr, taper_length, patches_out)

  :param blend_patch patches(\:):
  :param integer blend_mode:
  :param comm: Only present if compiled with MPI support, and of type `MPI_Comm` if using the F08 MPI interface or `integer` otherwise
  :param integer ierr: Will be zero on return if successful and non-zero otherwise
  :option integer taper_length:
  :option blend_patch patches_out(\:):

The parameters have the same meanings as in the C interface, with the addition of arrays of input and (optional) output patches of type :f:type:`blend_patch` that store patch-related values.

Values for the C enumerations (:c:enum:`AGDBlendMode` and :c:enum:`AGDTraceType`) are defined in the Fortran module.

.. f:subroutine:: deblend(patches, volumes, initial_factor, n_its, print_freq, comm, ierr, wavelets)

  :param deblend_patch patches(\:):
  :param volume volumes(\:):
  :param real initial_factor:
  :param integer n_its:
  :param integer print_freq:
  :param comm: Only present if compiled with MPI support, and of type `MPI_Comm` if using the F08 MPI interface or `integer` otherwise
  :param integer ierr: Will be zero on return if successful and non-zero otherwise
  :option wavelet wavelets(\:):

The deblending parameters are also similar to their C counterparts, again with the addition of arrays of types. This time patch-related values are stored in :f:type:`deblend_patch` arrays, while volume-related values use the :f:type:`volume` type and optional wavelets use the :f:type:`wavelet` type.

.. f:type:: blend_patch

  Patch for blending

  :f integer(kind=c_int) n_traces:
  :f integer(kind=c_int) n_times:
  :f integer(kind=c_long) shottimes(\:):
  :f integer(kind=c_int) channels(\:):
  :f integer(kind=c_int) trace_types(\:):
  :f real values(\:):

Even when using multiple space dimensions, the arrays that are included in :f:type:`blend_patch` and the other new types, should be one dimensional.

.. f:type:: deblend_patch

  Patch for deblending

  :f integer(kind=c_int) volume_idx:
  :f integer(kind=c_int) coords(\:):
  :f integer(kind=c_int) data_shape(\:):
  :f integer(kind=c_long) shottimes(\:):
  :f integer(kind=c_int) channels(\:):
  :f integer(kind=c_int) trace_types(\:):
  :f integer(kind=c_int) wavelet_idxs(\:):
  :f real values(\:): Data

As Fortran, like Julia, is a column-major language, some of the arrays in the Fortran interface, such as the `data_shape`, will also be provided in the opposite order to the C interface. As in the Julia interface, volume indexes should still be zero-indexed, however.

.. f:type:: volume

  :f integer(kind=c_int) n_dims:
  :f integer(kind=c_int) window_shape(\:):

.. f:type:: wavelet

  :f real values(\:): The source wavelet time series

C++
^^^

A C++ wrapper is not provided, as it should be possible to call the C interface directly. You will, however, have to apply `extern "C"` to the AGDeblend header file, as in `the C++ version of Example 1 <https://github.com/ar4/agdeblend/blob/main/examples/example_1.cpp>`_, if you compile the library with a C compiler. Some projects choose to include `extern "C"` in the header file itself so that this is not necessary, but AGDeblend does not do this so that the library can also be compiled with a C++ compiler, if desired.
