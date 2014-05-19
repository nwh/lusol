## LUSOL-Matlab interface

## Contents

General:

* `lusol_build.m`: script to generate prototype and thunk files
* `lusol.m`: function to obtain L and U factors from a factorization with LUSOL
* `lusol_obj.m`: interface to LUSOL for factorization, solves, and updates
* `lusol_test.m`: test suite for interface
* `test_data/`: data directory for test suite

Post build:

* `clusol.h`: header file for `clusol` library
* `libclusol_proto_glnxa64.m`: prototype file for Linux
* `libclusol_proto_maci61.m`: prototype file for Mac OS X
* `libclusol.so`: shared library for Linux
* `libclusol.dylib`: shared library for Mac OS X
* `libclusol_thunk_glnxa64.so`: thunk file for Linux
* `libclusol_thunk_maci64.dylib`: thunk file for Mac OS X
