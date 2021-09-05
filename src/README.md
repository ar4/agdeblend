# AGDeblend C Code

This directory contains the source code (written in C) for AGDeblend. The main file is `agdeblend.c`, which contains the public functions. A single translation unit approach is used, so this file uses `#include` to include the other source files during compilation, so this is the only file that needs to be passed to a compiler. The files ending with `_test.cpp` contain tests that use Google Test. The Makefile in the parents directory shows how to run them (you will first need to install Google Test).

See [the AGDeblend documentation](https://ausargeo.pages.dev/agdeblend) for more information.
