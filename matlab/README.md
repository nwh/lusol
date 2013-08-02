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

## Basic usage

In Matlab:

```
% get L and U factors
[L U P Q] = lusol(A);
```

See `help lusol`.

## Advanced usage

In Matlab:

```
% create lusol object
mylu = lusol_obj(A);

% solve with lusol object (ie x = A\b)
x = mylu.solveA(b); % x = A\b

% update factorization to replace a column
mylu.repcol(v,1);

% solve again with updated factorization
x1 = mylu.solveA(b);
```

See `help lusol_obj`.
