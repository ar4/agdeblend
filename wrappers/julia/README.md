# AGDeblend Julia Wrapper

This wrapper makes it easy to call AGDeblend from Julia.

To use it, please compile the [C AGDeblend](https://github.com/ar4/agdeblend) library first and make sure this library (or libraries, if you compile versions with different features) can be found by your operating system (such as by adding the directory to your `LD_LIBRARY_PATH` on Linux). You should then be able to use the library from Julia with:

```Julia
include(<path to agdeblend.jl>)
using .AGDeblend
```

If you wish to use MPI, please install the [MPI](https://www.juliapackages.com/p/mpi) package. You will probably need to configure it to use your system MPI library. This can be done with `julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'`. For more information, see the MPI package's [configuration documentation](https://juliaparallel.github.io/MPI.jl/stable/configuration/).

See [the AGDeblend documentation](https://ausargeo.pages.dev/agdeblend) for more information.
