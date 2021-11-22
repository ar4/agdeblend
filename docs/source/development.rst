Development
===========

This page provides an overview of how the AGDeblend code is structured. It is not necessary for users to read this, but might be useful for those who wish to understand or modify the code.

Operators
---------

Blending and deblending both involve linear operators. The blending operator is obviously important, and deblending, using the chosen approach of Fourier transformation, also requires a multidimensional Fourier transform operator (with windowing and padding), and, if source wavelets are provided, a wavelet convolution operator.

Blending
^^^^^^^^

The blending operator combines traces into a continuous record. AGDeblend provides versions that use the mean of overlapping samples or overwrites overlapping samples for working with recordings that are already blended, but the only version that is used within the deblending loop is the version that sums overlapping samples.

In theory, blending is quite a simple operation, as one simply needs to combine samples at the appropriate time and channel of the continuous record, but it is somewhat more complicated in AGDeblend in order to provide the disjoint continuous record, bad trace, and MPI features. Each channel of the continuous record is decomposed into intervals. Time samples are only stored for times within intervals. By having multiple intervals AGDeblend is able to support disjoint continuous records, as there can be an arbitrarily long gap between intervals.

To use the blending operator, we first need to set `blend_config`. This calculates the intervals that will be necessary to store the continuous record. It also calculates the samples of these intervals that are affected by traces labelled as bad and so should be muted in the continuous record. Furthermore, it determines the samples that overlap with intervals on other processes when using MPI. These overlapping samples will need to be shared between processes so that their value is the same on every process that has them.

In order to efficiently blend traces, we also need to set `blend_params` for each window of traces. This stores the index of the interval that each trace should be added to, as well as the number of samples into the interval to start from, and thus saves us from having to search through all of the intervals to find the appropriate one each time we want to apply the blending operator.

The adjoint of the blending operator copies samples from the continuous record into traces.

Fourier transform
^^^^^^^^^^^^^^^^^

In a typical dataset, the calculations performed by the Fourier transform operator will be the same for multiple windows. AGDeblend thus only creates a configuration, stored in `fk_config`, for each unique configuration, and then stores the index of the configuration to use for each window. The configuration contains FFTW plans and the taper values for each dimension.

The forward Fourier operator performs a complex-to-real multidimensional Fourier transform, to convert from the Fourier-domain model for the window into the time and space domain. A taper (described below) is applied, and a fraction of the length of the window in each dimension, corresponding to the region that is zero-padded in the adjoint operation, is discarded. Overlapping time windows are then summed to form traces. Each trace includes half a time window of extra samples at the beginning and end, compared to the corresponding input trace, that are tapered. This avoids hard edges in the continous record where traces start and stop, which would create noise over all frequencies in the Fourier transform.

The adjoint operator splits traces into time windows, adds the zero padding that was removed in the forward operator, applies the same taper, and performs a real-to-complex multidimensional Fourier transform. In order for this operator to be the adjoint of the forward operator, we then need to scale elements in the complex domain (except the first frequency component, and the last if the number of time samples is even) by two.

Wavelet convolution
^^^^^^^^^^^^^^^^^^^

In order to avoid having the effect of scaling components of the model, AGDeblend only uses the phase information from provided source wavelets. This is computed and saved in `source_configs` when the wavelet operator is being created. As multiple windows are likely to have the same shape, but may not necessarily use the same wavelets, the FFTW plan to transform a window in time that is needed to perform convolution with the wavelet, is stored separately. The wavelet operator parameters for each window, stored in `wavelet_params`, then specify which wavelet configurations and which FFTW plan to use for each window.

The forward and adjoint operators are identical except for a sign change. They Fourier transform the time dimension of the input window, multiply each trace by the phase of the wavelet that was specified for it, and then inverse Fourier transform the window.

Space Windows
-------------

AGDeblend splits patches into overlapping windows in every dimension when performing deblending. The time dimension is treated differently, however, as it is efficient to combine the time windows that make-up a trace before blending, so windows in the time dimension need to be processed together while windows in the other dimensions can be processed independently.

A list is constructed containing the space windows obtained from all of the patches in all of the volumes. Each iteration of deblending applies the forward and adjoint operators to each window in this list.

As the space windows overlap with each other, a taper is applied. The windows are tapered in every dimension except on edges that do not have any adjacent neighbouring windows. The taper is a Hann window taper, causing overlapping tapers to sum to one.

The space windows that a patch is decomposed into need to cover all of the traces in that patch. As the specified window shape might not divide evenly into the shape of the patch, some of the windows need to be bigger. The taper length must not be changed, however, so that the larger windows still overlap correctly with the smaller windows. Traces in the middle of the window that are not tapered are used to make windows larger.

Deblending Loop
---------------

The main work during deblending is performed in the loop that implements ISTA (iterative shrinkage and thresholding algorithm). This involves applying the forward operators to each window to calculate the continuously recorded, blended data that result from the current Fourier domain model. Calculating the residual between this and the recorded data (in continous record form), and applying the adjoint operators for each window will update the model. Finally, the shrink and threshold step is applied, which decreases the amplitude of each component of the Fourier domain model and zeros any values below the current threshold.

The threshold decays exponentially from the maximum initial value in the Fourier domain, multiplied by the initial factor. Exponential decay is used because we expect signal that passes the threshold to decrease the blending noise by a factor of the amplitude of that signal. The decay rate is chosen so that the threshold would reach `AGD_TARGET_LAMB` by the final iteration. In reality, however, the decay changes to a linear decrease to zero. The change happens at the iteration where the transition from the exponential to linear decrease is smooth.

AGD_TYPE
--------

The code is designed to support both single- and double-precision floating point values. To accomplish this, it uses several preprocessor macros. One of these is `AGD_TYPE`, which will be replaced by either `float` or `double` depending on which precision is specified during compilation. Others include `AGD_ONE`, replaced by `1.0f` or `1.0`, and `AGD_SQRT`, replace by `sqrtf` or `sqrt`.

Temporary Workspace
-------------------

Memory is needed to store intermediate values when applying the forward and adjoint operators to windows. This is provided by two pointers to allocated memory, called `temp_1` and `temp_2`. The memory that is allocated for these is chosen to be the largest size that will be needed for any of the windows. This approach means that these memory areas can be reused for every window, avoiding repeated costly memory allocations. When OpenMP is used, workspace is allocated for each thread so that they can work in parallel.

Source Files
------------

The source code files are in the `src <https://github.com/ar4/agdeblend/blob/main/src>`_ directory, with `agdeblend.c <https://github.com/ar4/agdeblend/blob/main/src/agdeblend.c>`_ containing the public functions. For this project, the advantages of a single translation unit were deemed to be greater than the disadvantages of this approach, so `agdeblend.c` uses `#include` to incorporate all of the other source files, and thus this is the only file that needs to be passed to a compiler.
