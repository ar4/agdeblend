# AGDeblend Examples

These examples demonstrate how to use AGDeblend for a variety of situations, calling it from C, Python, Julia, and Fortran (and C++ for Example 1).

Once you have compiled the AGDeblend code into a shared library, and installed the wrappers of the languages that require it, you can compile and run the examples as shown in the Makefile.

The Makefile recompiles the Fortran wrapper for each example, rather than compiling it once to a shared library. Both approaches should work equally well.
