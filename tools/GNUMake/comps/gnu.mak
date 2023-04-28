#
# Generic setup for using gcc
#
CXX = g++
CC  = gcc
FC  = gfortran
F90 = gfortran

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

gcc_version       := $(shell $(CXX) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;')
gcc_major_version := $(shell $(CXX) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;\..*;;')
gcc_minor_version := $(shell $(CXX) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

COMP_VERSION := $(gcc_version)

DEFINES += -DBL_GCC_VERSION='$(gcc_version)'
DEFINES += -DBL_GCC_MAJOR_VERSION=$(gcc_major_version)
DEFINES += -DBL_GCC_MINOR_VERSION=$(gcc_minor_version)

########################################################################

gcc_major_le_4 := $(shell expr $(gcc_major_version) \<= 4)
gcc_minor_lt_8 := $(shell expr $(gcc_minor_version) \< 8)
ifeq ($(gcc_major_le_4),1)
  ifeq ($(gcc_minor_lt_8),1)
    $(error GCC >= 4.8 required! Your version is $(gcc_version))
  endif
endif

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -Og -fno-inline -ggdb -Wshadow -Wall -Wno-sign-compare -ftrapv -Wno-unused-but-set-variable
  CFLAGS   += -Og -fno-inline -ggdb -Wshadow -Wall -Wno-sign-compare -ftrapv

  FFLAGS   += -Og -ggdb -pedantic -fcheck=all -fbacktrace -Wall -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv
  F90FLAGS += -Og -ggdb -pedantic -fcheck=all -fbacktrace -Wall -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv

  ifneq ($(gcc_major_version),$(filter $(gcc_major_version),4 5))
    CXXFLAGS += -Wnull-dereference
    CFLAGS += -Wnull-dereference
  endif

else

  CXXFLAGS += -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe -fopenmp
  CFLAGS   += -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe -fopenmp
  FFLAGS   += -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe -fopenmp
  F90FLAGS += -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe -fopenmp

endif


ifeq ($(USE_GPROF),TRUE)

  CXXFLAGS += -pg
  CFLAGS += -pg
  FFLAGS += -pg
  F90FLAGS += -pg

endif

########################################################################

ifeq ($(gcc_major_version),4)
  CXXFLAGS += -std=c++11
else ifeq ($(gcc_major_version),5)
  CXXFLAGS += -std=c++14
endif
CFLAGS     += -std=gnu99

FFLAGS   += -ffixed-line-length-none -fno-range-check -fno-second-underscore -J$(fmoddir) -I $(fmoddir)
F90FLAGS += -ffree-line-length-none -fno-range-check -fno-second-underscore -J$(fmoddir) -I $(fmoddir) -fimplicit-none

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(THREAD_SANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=thread
endif
ifeq ($(FSANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=address -fsanitize=undefined
endif

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -fopenmp
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

# ask gfortran the name of the library to link in.  First check for the
# static version.  If it returns only the name w/o a path, then it
# was not found.  In that case, ask for the shared-object version.
gfortran_liba  := $(shell $(F90) -print-file-name=libgfortran.a)
gfortran_libso := $(shell $(F90) -print-file-name=libgfortran.so)

ifneq ($(gfortran_liba),libgfortran.a)  # if found the full path is printed, thus `neq`.
  LIBRARY_LOCATIONS += $(dir $(gfortran_liba))
else
  LIBRARY_LOCATIONS += $(dir $(gfortran_libso))
endif

override XTRALIBS += -lgfortran -lquadmath
