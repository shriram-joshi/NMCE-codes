#-*- mode: makefile;-*-
CC=gcc
FC=gfortran
LD=gfortran

ifeq ($(loc),lap)
INC = #-I/home/vishrut.garg/tecplot/include -I/home/chris/usr/bin/openblas/include
LIB = #-L/home/vishrut.garg/tecplot/bin -L/home/chris/usr/bin/openblas/lib -lopenblas

else

INC = -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include
LIB = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
#INC = -I/opt/tecplot/current/include #-I/home/atom/a/anthonc/usr/bin/openblas/include
#LIB = -L/opt/tecplot/current/bin # -ltecio -lstdc++
#-L/home/atom/a/anthonc/usr/bin/openblas/lib -lopenblas

endif



ifeq ($(c),pg)
CC=pgcc
FC=pgf90

ifeq ($(opt),debug)
CFLAGS= -Mfree -C -g -Mchkfpstk -Mchkptr -traceback
FFLAGS= $(INC) -module Modules/ -Mfree -C -g -Mchkfpstk -Mchkptr -traceback
LDFLAGS= -module Modules/ -Mfree -C -g -Mchkfpstk -Mchkptr -traceback
else
#-O2 -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mvect=sse -Mcache_align -Mflushz -Mpre
-O2 -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mcache_align -Mflushz -Mpre
#CFLAGS= -Mfree -O2  -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mvect=sse -Mcache_align -Mflushz -Mpre
CFLAGS= -Mfree -O2  -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mcache_align -Mflushz -Mpre
#FFLAGS= -module Modules/ -Mfree -fast -Mvect=nosse
FFLAGS= -module Modules/ -Mfree -Mvect=nosse
#LDFLAGS= -module Modules/ -Mfree -O2 -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mvect=sse -Mcache_align -Mflushz -Mpre
LDFLAGS= -module Modules/ -Mfree -O2 -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mcache_align -Mflushz -Mpre

#CFLAGS= -Mfree -fast
#FFLAGS= -Mfree $(INC) -fast
#LDFLAGS= -Mfree -fast
endif

else ifeq ($(c),intel)
CC=icc
FC=ifort

ifeq ($(opt),debug)
#CFLAGS=-O3 -w -ipo -xhost
#FFLAGS= $(INC) -O3 -qopenmp -ipo -xhost -FR -u -module Modules -g -traceback -heap-arrays
#LDFLAGS= -O3 -qopenmp -ipo -xhost -module Modules -g -traceback -heap-arrays
CFLAGS=-O3 -w -xhost
#FFLAGS= $(INC) -O3 -qopenmp -xhost -FR -u -module Modules -g -traceback -heap-arrays
#LDFLAGS= -O3 -qopenmp -xhost -module Modules -g -traceback -heap-arrays
FFLAGS= $(INC) -O3 -qopenmp -xhost -FR -u -module Modules -V -O3 -debug -traceback -g -O2 -check bounds -C -heap-arrays
LDFLAGS= -O3 -qopenmp -xhost -module Modules -V -O3 -debug -traceback -g -O2 -check bounds -C -heap-arrays
#FFLAGS= $(INC) -O3 -qopenmp -xhost -FR -u -module Modules -g -check all -traceback -warn interfaces -fstack-protector -C -heap-arrays
#LDFLAGS= -O3 -qopenmp -xhost -module Modules -g -check all -traceback -warn interfaces -fstack-protector -C -heap-arrays
else
#CFLAGS=-O3 -w -ipo -xhost
#FFLAGS= $(INC) -O3 -pg -qopenmp -ipo -xhost -FR -u -module Modules -heap-arrays
#LDFLAGS= -O3 -pg -qopenmp -ipo -xhost -module Modules -heap-arrays
CFLAGS=-O3 -w -xhost
FFLAGS= $(INC) -O3 -qopenmp -xhost -FR -u -module Modules -heap-arrays
LDFLAGS= -O3 -qopenmp -xhost -module Modules -heap-arrays
endif

else
#CC=icc
#FC=ifort
#
#ifeq ($(opt),debug)
#CFLAGS=-O3 -w -ipo -xhost
#FFLAGS= $(INC) -O3 -qopenmp -ipo -xhost -FR -u -module Modules -g -traceback -heap_arrays
#LDFLAGS= -O3 -qopenmp -ipo -xhost -module Modules -g -traceback -heap_arrays
#else
#CFLAGS=-O3 -w -ipo -xhost
#FFLAGS= $(INC) -O3 -pg -qopenmp -ipo -xhost -FR -u -module Modules -heap_arrays
#LDFLAGS= -O3 -pg -qopenmp -ipo -xhost -module Modules -heap_arrays
#endif
CC=gcc
FC=gfortran
LD=gfortran

ifeq ($(opt),debug)

CFLAGS=-pg
FFLAGS = $(INC) -JModules -fcray-pointer -fimplicit-none -ffree-form -fbackslash -pg -g -fbounds-check -Wconversion -ffpe-trap=invalid,zero,overflow
LDFLAGS = -JModules -pg -fbounds-check -fopenmp -ffpe-trap=invalid,zero,overflow

else ifeq ($(opt),debug2)

CFLAGS=-O2 -fopenmp
FFLAGS=$(INC) -O3 -g -JModules -ffree-form -fopenmp -march=native -funroll-loops -ffpe-trap=invalid,zero,overflow  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-arrays
LDFLAGS= -O3 -g -JModules -fopenmp -ffree-form -march=native -funroll-loops -ffpe-trap=invalid,zero,overflow  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-array

else ifeq ($(opt),prof)

CFLAGS=-O2 -fopenmp
FFLAGS=$(INC) -O3 -pg -JModules -ffree-form -fopenmp -march=native -funroll-loops  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-arrays
LDFLAGS= -O3 -pg -JModules -fopenmp -ffree-form -march=native -funroll-loops  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-array

else ifeq ($(opt),fast)

CFLAGS=-O2 -fopenmp
FFLAGS=$(INC) -O3 -JModules -ffree-form -fopenmp -march=native -funroll-loops -ffast-math -fno-protect-parens -flto
LDFLAGS= -O3 -JModules -fopenmp -ffree-form -march=native -funroll-loops -ffast-math -fno-protect-parens -flto

else ifeq ($(opt),asan)

CFLAGS=-O2 -fopenmp -fsanitize=address
FFLAGS=$(INC) -O3 -JModules -ffree-form -fopenmp -march=native -funroll-loops -fsanitize=address
LDFLAGS= -O3 -JModules -fopenmp -ffree-form -march=native -funroll-loops -fsanitize=address

else

CFLAGS=-O2 -fopenmp
FFLAGS=$(INC) -O3 -JModules -ffree-form -fopenmp -march=native -funroll-loops
LDFLAGS= -O3 -JModules -fopenmp -ffree-form -march=native -funroll-loops

endif

endif

SRCDIR=Sources
OBJDIR=Objects

#Add sources as $(SRCDIR)/file.f90 or $(SRCDIR)/file.c
CSOURCES=
FSOURCES= $(SRCDIR)/set_kd.f90 $(SRCDIR)/io.f90 $(SRCDIR)/data.f90 $(SRCDIR)/library.f90 $(SRCDIR)/mesh.f90 $(SRCDIR)/assembly.f90  $(SRCDIR)/Bandedsolver.f90 $(SRCDIR)/main.f90 ##########################################

FOBJECTS= $(FSOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
COBJECTS= $(CSOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
COBJECTSCMD=$(CSOURCESCMD:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

EXECUTABLE= output

all: $(CSOURCES) $(FSOURCES) $(EXECUTABLE)

cmd: $(FSOURCES) $(CSOURCESCMD) $(EXECUTABLECMD)

$(EXECUTABLE): $(FOBJECTS) $(COBJECTS)
	$(FC) $(LDFLAGS) $(FOBJECTS) $(COBJECTS) -o $@ $(LIB)

$(FOBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(COBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm Modules/*mod $(EXECUTABLE) $(FOBJECTS) $(COBJECTS)

wipe:
	rm Modules/*mod $(EXECUTABLE) $(FOBJECTS) $(COBJECTS)
