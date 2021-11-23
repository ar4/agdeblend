from cffi import FFI
from agdeblend_libs import lib_dir, lib_exists

cdef = """
enum AGDTraceType { AGDLive, AGDBad, AGDMissing };
enum AGDBlendMode { AGDBlendSum, AGDBlendMean, AGDBlendOverwrite };
int agd_deblend(int n_patches, int const *volumes, int const *n_dims,
                int const *const *window_shapes, int const *const *coords,
                int const *const *data_shapes,
                long int const *const *shottimes,
                int const *const *channels,
                enum AGDTraceType const *const *trace_types,
                int const *wavelet_lengths,
                int const *const *wavelet_idxs,
                float const *const *wavelets,
                float initial_factor, int n_its, int print_freq,
                MPI_Comm comm,
                float *const *data);
int agd_blend(int n_patches, int const *n_traces,
              int const *n_times,
              long int const *const *shottimes,
              int const *const *channels,
              enum AGDTraceType const *const *trace_types,
              float *const *data, enum AGDBlendMode blend_mode,
              int taper_length, int n_patches_out,
              int const *n_traces_out,
              int const *n_times_out,
              long int const *const *shottimes_out,
              int const *const *channels_out,
              enum AGDTraceType const *const *trace_types_out,
              MPI_Comm comm,
              float *const *data_out);
"""

if lib_exists["agdeblend"]:
    ffibuilder = FFI()
    lib_cdef = cdef.replace("MPI_Comm comm,", "")
    ffibuilder.cdef(lib_cdef)
    ffibuilder.set_source(
        "agdeblend._agdeblend_ffi",
        lib_cdef,
        library_dirs=[lib_dir],
        libraries=["agdeblend"],
    )

if lib_exists["agdeblend_double"]:
    ffibuilder_double = FFI()
    lib_cdef = cdef.replace("MPI_Comm comm,", "").replace("float", "double")
    ffibuilder_double.cdef(lib_cdef)
    ffibuilder_double.set_source(
        "agdeblend._agdeblend_double_ffi",
        lib_cdef,
        library_dirs=[lib_dir],
        libraries=["agdeblend_double"],
    )


def define_mpi_comm(ffi):
    from mpi4py import MPI

    if MPI._sizeof(MPI.Comm) == ffi.sizeof("int"):
        ffi.cdef("typedef int MPI_Comm;")
    elif MPI._sizeof(MPI.Comm) == ffi.sizeof("void *"):
        ffi.cdef("typedef void * MPI_Comm;")
    else:
        raise RuntimeError("MPI_Comm is of unexpected type")


if lib_exists["agdeblend_mpi"]:
    ffibuilder_mpi = FFI()
    define_mpi_comm(ffibuilder_mpi)
    ffibuilder_mpi.cdef(cdef)
    ffibuilder_mpi.set_source(
        "agdeblend._agdeblend_mpi_ffi",
        "#include <mpi.h>" + cdef,
        library_dirs=[lib_dir],
        libraries=["agdeblend_mpi"],
    )

if lib_exists["agdeblend_mpi_double"]:
    ffibuilder_mpi_double = FFI()
    define_mpi_comm(ffibuilder_mpi_double)
    lib_cdef = cdef.replace("float", "double")
    ffibuilder_mpi_double.cdef(lib_cdef)
    ffibuilder_mpi_double.set_source(
        "agdeblend._agdeblend_mpi_double_ffi",
        "#include <mpi.h>" + lib_cdef,
        library_dirs=[lib_dir],
        libraries=["agdeblend_mpi_double"],
    )

if __name__ == "__main__":
    if lib_exists["agdeblend"]:
        ffibuilder.compile(verbose=True)
    if lib_exists["agdeblend_double"]:
        ffibuilder_double.compile(verbose=True)
    if lib_exists["agdeblend_mpi"]:
        ffibuilder_mpi.compile(verbose=True)
    if lib_exists["agdeblend_mpi_double"]:
        ffibuilder_mpi_double.compile(verbose=True)
