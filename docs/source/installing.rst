Installing
==========

AGDeblend is written in C and so needs to be compiled. The only two requirements are a C compiler and a Fourier transform library with an FFTW3 interface. If you wish to use the optional MPI and OpenMP features for parallel processing, you will also need libraries that provide those. Once you have these, compiling and using AGDeblend should be easy.

After AGDeblend has been installed, you might want to use the provided wrappers to call it from other languages. Instructions are below.

AGDeblend Libraries
-------------------

You will first need to download the source code. You can do this by downloading the `latest release <https://github.com/ar4/agdeblend/releases/latest>`_ and extracting it.

Only one source file needs to be compiled: ``src/agdeblend.c``. This will automatically include the other source files. The recommended approach is to compile AGDeblend into a shared library that you can then link with your own code. For a single-threaded installation using single-precision floating point values and the FFTW3 library, a typical compilation command when using the GCC compiler is::

  mkdir lib
  gcc -shared -fPIC -O2 src/agdeblend.c -Iinclude -lfftw3f -lm -o lib/libagdeblend.so

You are recommended to either copy the output library to a directory that your operating will search for libraries, or add the ``lib`` directory that the library was saved in to the list of directories that will be searched. On Linux the latter can be done by adding the path of this directory to the ``LD_LIBRARY_PATH`` environment variable (do it in your startup script so that it will always be there).

A higher optimisation level, such as `O3` or `Ofast`, might result in better performance on your system. If you know which processors the library will execute on, enabling features supported by all of those processors, perhaps including AVX2, might also produce a performance improvement. If you will run on the same processor type as used during compilation, you can enable the features that it supports using `-march=native`.

If you want to use double-precision instead of single-precision, define the ``AGD_DOUBLE`` macro when compiling, which can be done in most compilers by adding ``-DAGD_DOUBLE`` to the compilation command. You will also need to link with the double-precision FFTW3 library instead of the single-precision one, which can be done by replacing ``-lfftw3f`` with ``-lfftw3``.

To use more than one processing core, you can enable OpenMP. This may cause deblending to run somewhat faster, but you are unlikely to see a linear speedup. It is, however, easy to use, as you simply need to define the ``AGD_THREADS`` macro when compiling and enable your compiler's OpenMP support. You do not need to make any changes to how you use AGDeblend. With most compilers, if you have OpenMP installed, this can be done by adding ``-DAGD_THREADS -fopenmp`` to the compilation command.

Enabling MPI support, which allows you to run blending and deblending on multiple processes, on one or more processors, is also simple if you already have MPI installed. You just need to define the ``AGD_MPI`` macro and link with your MPI library. The latter is usually most easily accomplished by using the ``mpicc`` wrapper that is typically included with MPI installations, in place of your regular compiler.

AGDeblend should work with any Fourier transform library that supports the FFTW3 interface. In addition to FFTW3 itself, another popular choice is Intel's MKL. To do this you will replace ``-lfftw3f`` in the above command. What to replace it with depends on your particular setup. The easiest way to determine this is probably to use Intel's online Link Line Advisor. The sequential (non-threaded) version of MKL is recommended for use with AGDeblend.

A Makefile is provided in the base directory of AGDeblend. For a typical setup using GCC and FFTW3, this should usually be sufficient.

AGDeblend is written so that it should also be possible to compile it using a C++ compiler, if necessary.

For your own code that uses AGDeblend, simply link with the installed library, for example with ``-lagdeblend``.

Wrappers
--------

Once you have the AGDeblend library installed, and have added it to the paths that you system checks for libraries, you may call it through the provided Python, Julia, or Fortran wrappers if you prefer.

Python
^^^^^^

To install the Python wrapper, enter the ``wrappers/python`` directory, and run ``AGD_LIB_DIR=<absolute path to installed AGDeblend library> pip install .``. This assumes that the single-precision non-MPI library is called ``agdeblend`` (with or without ``lib`` at the front), the double-precision version adds ``_double`` afterwards, and the MPI version(s) add ``_mpi`` at the end. If you did not install all of these versions of the library, the Python wrapper support for those features will not be installed. Once this has completed, you should be able to ``import agdeblend`` into your own code.

Julia
^^^^^

The Julia wrapper does not currently require installation.

Fortran
^^^^^^^

The Fortran wrapper consists of a Fortran file (using the Fortran 2003 standard) and a C file.

You should compile the C file first. The compilation is similar to that of the AGDeblend library, with the same ``AGD_DOUBLE`` and ``AGD_MPI`` macros, but ``AGD_THREADS`` will have no effect and linking with FFTW and OpenMP is not required. As with AGDeblend, you will need to include the path to the ``agdeblend.h`` header file, which is in the ``include`` directory inside AGDeblend's base directory. You should only compile this file into a shared object, which can be achieved with most compilers using ``-c -o agdeblend_c.o``.

You are now ready to compile the Fortran component of the wrapper. Compilation is the same as for the C file, but with the addition of an extra macro ``AGD_F08`` which you should define (in addition to ``AGD_MPI``) if you are using the Fortran 2008 MPI interface. You will also need to add the shared object produced by compiling the C file to the list of input files, and to link with the AGDeblend library. You can create a shared library as the result (with ``-shared -fPIC -o libagdeblend_f.so``, for example). As with the original AGDeblend library, you should make sure that your operating system can find this library, either by putting it in a standard location or by adding the path to the list of directories that are searched. When compiling your own Fortran code, you should then link with this new library (``-lagdeblend_f``) and you will also need to tell your compiler where to find the module file that was created (with GCC, this is done with ``-J <path to module file>``).
