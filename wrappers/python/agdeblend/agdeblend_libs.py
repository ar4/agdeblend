import os

try:
    lib_dir = os.path.abspath(os.environ["AGD_LIB_DIR"])
except KeyError:
    raise RuntimeError(
        "Please set the AGD_LIB_DIR environment variable to the path "
        "of the directory containing the AGDeblend library files"
    )

lib_exists = {}
n_found = 0
for libname in [
    "agdeblend",
    "agdeblend_double",
    "agdeblend_mpi",
    "agdeblend_mpi_double",
]:
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
            found = True
            n_found += 1
            break
    lib_exists[libname] = found
if n_found == 0:
    raise RuntimeError("No AGDeblend compiled library was found")
