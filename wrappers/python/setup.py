import os
import setuptools
from setuptools import setup
from agdeblend_libs import lib_exists

cffi_modules = []
if lib_exists["agdeblend"]:
    cffi_modules.append("agdeblend_build.py" + ":ffibuilder")
if lib_exists["agdeblend_double"]:
    cffi_modules.append("agdeblend_build.py" + ":ffibuilder_double")
if lib_exists["agdeblend_mpi"]:
    cffi_modules.append("agdeblend_build.py" + ":ffibuilder_mpi")
if lib_exists["agdeblend_mpi_double"]:
    cffi_modules.append("agdeblend_build.py" + ":ffibuilder_mpi_double")

setup(
    name="agdeblend",
    version="1.0.1",
    description="Seismic data blending and deblending",
    url="https://github.com/ar4/agdeblend",
    author="Alan Richardson",
    author_email="alan@ausargeo.com",
    license="GPLv3",
    project_urls={
      "Documentation": "https://ausargeo.pages.dev/agdeblend"
    },
    packages=setuptools.find_packages(),
    setup_requires=["cffi>=1.0.0"],
    extras_require={"mpi": ["mpi4py"]},
    cffi_modules=cffi_modules,
    install_requires=["cffi>=1.0.0", "numpy"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
