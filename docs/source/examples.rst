Examples
========

The easiest way to learn how to use AGDeblend is to see examples, so several are provided for a range of common situations, in each of the supported programming languages.

In each case, the code assumes that the input arrays have already been created `make_data.c <https://github/com/ar4/agdeblend/blob/master/examples/make_data.c>`_ can be used), loads them, performs the task, and then saves the output to disk.

The code for all of the examples can be found in the `examples <https://github/com/ar4/agdeblend/blob/master/examples/>`_ directory.

.. toctree::

  Simple synthetic blending of unblended data <example_1>
  Adjust already blended data to have longer traces <example_2>
  Blending using MPI <example_3>
  Simple deblending of one patch <example_4>
  Deblending of two patches <example_5>
  Deblending with MPI <example_6>
  Deblending with wavelet convolution <example_7>
  Deblending two volumes <example_8>
  Deblending an irregularly shaped volume <example_9>
