# AGDeblend Fortran Wrapper

This wrapper makes it easy to call AGDeblend from Fortran.

To use it, please compile the [C AGDeblend](https://github.com/ar4/agdeblend) library first and make sure this library (or libraries, if you compile versions with different features) can be found by your operating system (such as by adding the directory to your `LD_LIBRARY_PATH` on Linux). You should then compile `agdeblend_c.c` to a shared object, and then compile this object file with `agdeblend.F90` to a shared library.

See [the AGDeblend documentation](https://ausargeo.pages.dev/agdeblend) for more information.
