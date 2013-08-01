## LUSOL

LUSOL maintains LU factors of a square or rectangular sparse matrix.

This repository provides LUSOL source code and a Matlab interface.

The website for LUSOL is: [http://www.stanford.edu/group/SOL/software/lusol.html][LUSOL]

The code is distributed under the terms of the [Common Public License][CPL].

  [LUSOL]: http://www.stanford.edu/group/SOL/software/lusol.html
  [CPL]: http://www.ibm.com/developerworks/library/os-cpl.html

## Contents

* `gen/`: code generation scripts and specification files
* `matlab/`: Matlab interface code
* `src/`: LUSOL fortran code
* `LICENSE`: [Common Public License][CPL]
* `README.md`: this file
* `makefile`: GNU Make file to build interface

## Matlab install

Add `lusol/matlab` to your Matlab path.  This may be done with Matlab's [`addpath`][ADDPATH] function.

If the interface has not been built, please follow the directions below.

  [ADDPATH]: http://www.mathworks.com/help/matlab/ref/addpath.html

## Build

From the `lusol` directory:

```
$ make
```

Linux output from 2013-07-31:

```
$ make
gcc  -fPIC -c src/clusol.c -o src/clusol.o
gfortran -fPIC -Jsrc -c src/lusol_precision.f90 -o src/lusol_precision.o
gfortran -fPIC -Jsrc -c src/lusol.f90 -o src/lusol.o
gfortran -fPIC -fdefault-integer-8 -c src/lusol_util.f -o src/lusol_util.o
gfortran -fPIC -fdefault-integer-8 -c src/lusol6b.f -o src/lusol6b.o
gfortran -fPIC -fdefault-integer-8 -c src/lusol7b.f -o src/lusol7b.o
gfortran -fPIC -fdefault-integer-8 -c src/lusol8b.f -o src/lusol8b.o
gcc src/clusol.o src/lusol.o src/lusol_precision.o src/lusol_util.o src/lusol6b.o src/lusol7b.o src/lusol8b.o -o src/libclusol.so -shared -lgfortran
mkdir -p lib
cp src/libclusol.so lib/libclusol.so
mkdir -p include
cp src/clusol.h include/clusol.h
cp src/libclusol.so src/clusol.h ./matlab/
matlab -nojvm -nodisplay -r "cd matlab; lusol_build; exit"

                            < M A T L A B (R) >
                  Copyright 1984-2011 The MathWorks, Inc.
                    R2011b (7.13.0.564) 64-bit (glnxa64)
                              August 13, 2011


To get started, type one of these: helpwin, helpdesk, or demo.
For product information, visit www.mathworks.com.
```

## Test

From the lusol directory:

```
$ make matlab_test
```

Linux output from 2013-07-31:

```
$ make matlab_test
matlab -nojvm -nodisplay -r "cd matlab; lusol_test; exit"

                            < M A T L A B (R) >
                  Copyright 1984-2011 The MathWorks, Inc.
                    R2011b (7.13.0.564) 64-bit (glnxa64)
                              August 13, 2011


To get started, type one of these: helpwin, helpdesk, or demo.
For product information, visit www.mathworks.com.

test_addcol_lhr02: passed
test_addrow_lhr02: passed
test_delcol_lhr02: passed
test_delrow_lhr02: passed
test_factorize_lhr02: passed
test_mulA_lhr02: passed
test_mulAt_lhr02: passed
test_r1mod_lhr02: passed
test_repcol_lhr02: passed
test_reprow_lhr02: passed
test_solveA_lhr02: passed
test_solveAt_lhr02: passed
```

## Authors

* **LUSOL Fortran code**: [Michael Saunders][MS]
* **Matlab interface**: [Nick Henderson][NWH]

  [MS]: http://www.stanford.edu/~saunders/
  [NWH]: http://www.stanford.edu/~nwh/
