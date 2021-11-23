import os

try:
    from mpi4py import MPI

    mpi4py_installed = True
except ModuleNotFoundError:
    mpi4py_installed = False

try:
    lib_dir = os.path.abspath(os.environ["AGD_LIB_DIR"])
except KeyError:
    raise RuntimeError(
        "Please set the AGD_LIB_DIR environment variable to the path "
        "of the directory containing the AGDeblend library files"
    )

lib_exists = {}
n_found = 0
for lib in [
    {"libname": "agdeblend", "mpi": False},
    {"libname": "agdeblend_double", "mpi": False},
    {"libname": "agdeblend_mpi", "mpi": True},
    {"libname": "agdeblend_mpi_double", "mpi": True},
]:
    libname = lib["libname"]
    possible_filenames = [
        "lib" + libname + ".so",
        libname + ".so",
        "lib" + libname + ".dll",
        libname + ".dll",
        "lib" + libname + ".dylib",
        libname + ".dylib",
    ]
    found = False
    for possible_filename in possible_filenames:
        if os.path.exists(os.path.join(lib_dir, possible_filename)):
            if lib["mpi"] and not mpi4py_installed:
                print("Please install mpi4py before agdeblend to use MPI")
                break
            else:
                found = True
                n_found += 1
                break
    lib_exists[libname] = found
if n_found == 0:
    raise RuntimeError("No AGDeblend compiled library was found")
