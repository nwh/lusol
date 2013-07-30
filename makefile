# Remove default rules
.SUFFIXES:

# detect operating system
OSLOWER := $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])
DARWIN := $(strip $(findstring darwin, $(OSLOWER)))

# C compiler
CC := gcc
CPPFLAGS :=
CFLAGS := -fPIC

# Fortran 90 compiler
F90C := gfortran
F90FLAGS := -fPIC -Jsrc

# Fortran 77 compiler
F77C := gfortran
F77FLAGS := -fPIC -fdefault-integer-8

# Matlab
ML := matlab
MLFLAGS := -nojvm -nodisplay
ifneq ($(DARWIN),)
  # settings for mac os x
  MLARCH := maci64
else
  # settins for linux
  MLARCH := glnxa64
endif

# Linker
LD := gcc
ifneq ($(DARWIN),)
  # settings for mac os x
  LIB_SUFFIX := dylib
  LDFLAGS := -dynamiclib
  LDFLAGS += -L/Applications/MATLAB_R2011b.app/bin/maci64
  LDFLAGS += -L/usr/local/Cellar/gfortran/4.8.1/gfortran/lib/gcc/x86_64-apple-darwin12.3.0/4.8.1/../../../
  LDFLAGS += -lmwblas -lgfortran
else
  # settins for linux
  LIB_SUFFIX := so
  LDFLAGS := -shared -lgfortran
endif

# list of files required by matlab
MATLAB_FILES := \
  matlab/libclusol.$(LIB_SUFFIX) \
  matlab/clusol.h \
  matlab/libclusol_proto_$(MLARCH).m \
  matlab/libclusol_thunk_$(MLARCH).$(LIB_SUFFIX)

# list of interface specification files
INTERFACE_FILES := \
  gen/interface.py \
  gen/interface_files.org \
  gen/lu1fac.org \
  gen/lu6mul.org \
  gen/lu6sol.org \
  gen/lu8adc.org \
  gen/lu8adr.org \
  gen/lu8dlc.org \
  gen/lu8dlr.org \
  gen/lu8mod.org \
  gen/lu8rpc.org \
  gen/lu8rpr.org

# list of F77 code files
F77_FILES := \
	src/lusol_util.f \
  src/lusol6b.f \
  src/lusol7b.f \
  src/lusol8b.f

F77_OBJ := $(patsubst %.f,%.o,$(filter %.f,$(F77_FILES)))

# default target to build everything
.PHONY: all
all: lib/libclusol.$(LIB_SUFFIX) include/clusol.h $(MATLAB_FILES)

# pattern to compile fortran 77 files
$(F77_OBJ) : %.o : %.f
	$(F77C) $(F77FLAGS) -c $< -o $@

src/lusol_precision.o: src/lusol_precision.f90
	$(F90C) $(F90FLAGS) -c $< -o $@

src/lusol_precision.mod: src/lusol_precision.o

src/lusol.o: src/lusol.f90 src/lusol_precision.mod
	$(F90C) $(F90FLAGS) -c $< -o $@

src/lusol.mod: src/lusol.o

src/clusol.h: $(INTERFACE_FILES)
	./gen/interface.py -i gen/interface_files.org -o $@ -t header

src/clusol.c: $(INTERFACE_FILES)
	./gen/interface.py -i gen/interface_files.org -o $@ -t source

src/clusol.o: src/clusol.c src/clusol.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

src/libclusol.$(LIB_SUFFIX): src/clusol.o src/lusol.o src/lusol_precision.o $(F77_OBJ)
	$(LD) $^ -o $@ $(LDFLAGS)

include/clusol.h: src/clusol.h
	mkdir -p include
	cp $< $@

lib/libclusol.$(LIB_SUFFIX): src/libclusol.$(LIB_SUFFIX)
	mkdir -p lib
	cp $< $@

$(MATLAB_FILES): src/libclusol.$(LIB_SUFFIX) src/clusol.h
	cp src/libclusol.$(LIB_SUFFIX) src/clusol.h ./matlab/
	$(ML) $(MLFLAGS) -r "cd matlab; lusol_build; exit"

.PHONY: matlab
matlab: $(MATLAB_FILES)

.PHONY: matlab_test
matlab_test: $(MATLAB_FILES)
	$(ML) $(MLFLAGS) -r "cd matlab; lusol_test; exit"

.PHONY: clean
clean:
	$(RM) src/*.o
	$(RM) src/*.$(LIB_SUFFIX)
	$(RM) src/*.mod
	$(RM) lib/libclusol.so
	$(RM) lib/libclusol.$(LIB_SUFFIX)
	$(RM) include/clusol.h
	$(RM) $(MATLAB_FILES)

.PHONY: clean_gen
clean_gen:
	$(RM) src/clusol.h
	$(RM) src/clusol.c

# print helper
#print-%:
#	@echo $* := $($*)
