## LUSOL fortran code

This directory contains LUSOL fortran 77 and fortran 90 code.

## Contents

* `clusol.c`: generated interface code
* `clusol.h`: generated interface header
* `lusol6b.f`: contains lu6mul
* `lusol7b.f`: contains lusol helper functions
* `lusol8b.f`: contains lusol update functions (not lu8rpc)
* `lusol.f90`: fortran 90 code for factorization, solve, and replace column
* `lusol_precision.f90`: fortran 90 module to specify precision
* `lusol.txt`: some LUSOL documentation and history
* `lusol_util.f`: fortran 77 code for factorization, solve, and replace column
* `updates/`: (2016-01-23) a directory containing updates to LUSOL, will be
  factored in to the Matlab interface when time permits.
