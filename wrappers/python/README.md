# AGDeblend Python Wrapper

This wrapper makes it easy to call AGDeblend from Python.

To use it, please compile the [C AGDeblend](https://github.com/ar4/agdeblend) library first, make sure this library (or libraries, if you compile versions with different features) can be found by your operating system (such as by adding the directory to your `LD_LIBRARY_PATH` on Linux), and then run `AGD_LIB_DIR=<path to directory with compiled AGDeblend library> pip install agdeblend`. You should then be able to run `import agdeblend` to use AGDeblend from Python.

If you wish to use MPI, please install [mpi4py](https://github.com/mpi4py/mpi4py) before agdeblend. This can be done automatically by adding `[mpi]` to the installation command. You will probably also need to specify that `mpicc` should be used as the compiler, so the installation command becomes `AGD_LIB_DIR=<path to AGDeblend libraries> CC=mpicc pip install agdeblend[mpi]`.

See [the AGDeblend documentation](https://ausargeo.pages.dev/agdeblend) for more information.
