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
* `src/`: LUSOL Fortran code
* `LICENSE`: [Common Public License][CPL]
* `makefile`: GNU Make file to build interface

## Download and install

Pre-built downloads are available on the Github [release][RELEASE] page.

Installation simply requires adding the `matlab` subdirectory to your Matlab
path.  This may be done with Matlab's [`addpath`][ADDPATH] function.

If the interface has not been built, please follow the directions below.

  [RELEASE]: https://github.com/nwh/lusol/releases
  [ADDPATH]: http://www.mathworks.com/help/matlab/ref/addpath.html

## Build

### Environment

The build environments as of 2014-05-18 are:

- Fedora 21 & Matlab 2013b
- Mac OS X 10.9 & Matlab 2014a

Building the LUSOL interface in other environments may require modification of
`makefile` and `matlab/lusol_build.m`.

### Requirements

Linux:

* `make`
* `gcc`
* `gfortran`
* Matlab

Mac:

* [Xcode][XC] for `clang` and `make`
* `gfortran-4.3` (possibly via [Homebrew][HB])
* Matlab

  [HB]: http://brew.sh/
  [XC]: http://itunes.apple.com/us/app/xcode/id497799835

Notes:

* The `matlab` binary must be on the system `PATH`.
* `python3` is required to generate the interface code.  However, the interface
  code is pre-generated and included in the repository.
* It may be necessary to launch Xcode and accept the license agreement before
  building the interface.

### Setup `mex`

Matlab's `mex` compiler driver must be configured to use the appropriate `C`
compiler.  This can be achieved by executing `mex -setup` from the Matlab prompt
or operating system terminal.  On Linux the selected compiler must be the
correct version of `gcc`.  On Mac OS X 10.9 the selected compiler must be
`clang`.  It is a good idea to match the compiler suggested on the Mathworks
[supported compilers][MC] page.

  [MC]: http://www.mathworks.com/support/compilers/

### Steps

From the base directory:

```
$ make
```

Mac OS X output from 2014-05-18:

```
$ make
clang  -fPIC -c src/clusol.c -o src/clusol.o
gfortran-4.3 -fPIC -Jsrc -O3 -c src/lusol_precision.f90 -o src/lusol_precision.o
gfortran-4.3 -fPIC -Jsrc -O3 -c src/lusol.f90 -o src/lusol.o
gfortran-4.3 -fPIC -fdefault-integer-8 -O3 -c src/lusol_util.f -o src/lusol_util.o
gfortran-4.3 -fPIC -fdefault-integer-8 -O3 -c src/lusol6b.f -o src/lusol6b.o
gfortran-4.3 -fPIC -fdefault-integer-8 -O3 -c src/lusol7b.f -o src/lusol7b.o
gfortran-4.3 -fPIC -fdefault-integer-8 -O3 -c src/lusol8b.f -o src/lusol8b.o
clang src/clusol.o src/lusol.o src/lusol_precision.o src/lusol_util.o src/lusol6b.o src/lusol7b.o src/lusol8b.o -o src/libclusol.dylib -dynamiclib -L/Applications/MATLAB_R2014a.app/bin/maci64 -lmwblas -L/Applications/MATLAB_R2014a.app/sys/os/maci64 -lgfortran
mkdir -p lib
cp src/libclusol.dylib lib/libclusol.dylib
mkdir -p include
cp src/clusol.h include/clusol.h
cp src/libclusol.dylib src/clusol.h ./matlab/
matlab -nojvm -nodisplay -r "cd matlab; lusol_build; exit"

                            < M A T L A B (R) >
                  Copyright 1984-2014 The MathWorks, Inc.
                     R2014a (8.3.0.532) 64-bit (maci64)
                             February 11, 2014


To get started, type one of these: helpwin, helpdesk, or demo.
For product information, visit www.mathworks.com.
```

### Test

From the base directory:

```
$ make matlab_test
```

Mac OS X output from 2014-05-18:

```
$ make matlab_test
matlab -nojvm -nodisplay -r "cd matlab; lusol_test; exit"

                            < M A T L A B (R) >
                  Copyright 1984-2014 The MathWorks, Inc.
                     R2014a (8.3.0.532) 64-bit (maci64)
                             February 11, 2014


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

### Notes

The basic requirements to build LUSOL are GNU `make`, `gfortran`, a C compiler,
and Matlab.  The build works with `gcc` on Linux and `clang` on Mac OS X.  It
may be necessary to modify the compiler variables in the `makefile` (`CC`,
`F90C`, and `F77C`) depending on the operating system and environment.

The `matlab` executable must be on the system path.  On Mac OS X with Matlab
R2014a this is achieved with:

```
$ export PATH=/Applications/MATLAB_R2014a.app/bin:$PATH
```

The `makefile` may have to be modified on Mac OS X depending on versions of
Matlab and `gfortran`.  The [`LDFLAGS`][LDFLAGS] need to point to the proper
directories for Matlab `blas` and `gfortran` shared libraries.

  [LDFLAGS]: https://github.com/nwh/lusol/blob/master/makefile#L58

## Authors

* **LUSOL Fortran code**: [Michael Saunders][MS]
* **Matlab interface**: [Nick Henderson][NWH]

  [MS]: http://www.stanford.edu/~saunders/
  [NWH]: http://www.stanford.edu/~nwh/
