#
# makefile
#

BINROOT = $(CURDIR)
FC=undefined
CC=undefined
FFLAGS=undefined
LDFLAGS=undefined
TYPE=undefined
EX=main.x
EXP=main_mpi.x
EXE = $(BINROOT)/$(EX)
EXEP = $(BINROOT)/$(EXP)
SHELL=/bin/sh

def:	all

all:
	@echo "Error - you should specify a target machine!"
	@echo "Possible choices are:              "
	@echo "                                   "
	@echo "nvfortran                          "
	@echo "nvfortran-mpi                      "
	@echo "gnu-debug                          "
	@echo "gnu-mpi-debug                      "
	@echo "gnu-real                           "
	@echo "gnu-mpi-real                       "
	@echo "gnu-double                         "
	@echo "gnu-mpi-double                     "
	@echo "                                   "
	@echo "Possible choices for debugging are:"
	@echo "                                   "
	@echo "                                   "
	@echo "Please examine Makefile for further details "

help: all
	
nvfortran:
	$(MAKE) CC=nvcc \
	FC=nvfortran \
	CFLAGS="-O3 -I $(CURDIR) -c" \
	FFLAGS="-O3 -cpp -acc -gpu=cc75 -DPRC=4 -DLATTICE=D3Q27 -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-O3 -cpp -acc -gpu=cc75 -DPRC=4 -DVERBOSE -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

nvfortran-mpi:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O1 -I $(CURDIR) -c" \
	FFLAGS="-O1 -cpp -acc -gpu=cc75 -DPRC=4 -DLATTICE=D3Q27 -DMPI -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-O1 -cpp -acc -gpu=cc75 -DPRC=4 -DMPI -DVERBOSE -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi

gnu-debug:
	$(MAKE) CC=gcc \
	FC=gfortran \
	CFLAGS="-O1 -g -I $(CURDIR) -c" \
	FFLAGS="-O1 -g -fcheck=all -cpp -DPRC=4 -DLATTICE=D3Q27 -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-O1 -g -fcheck=all -cpp -DPRC=4 -DVERBOSE -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

gnu-mpi-debug:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O1 -g -I $(CURDIR) -c" \
	FFLAGS="-O1 -g -fcheck=all -cpp -DPRC=4 -DLATTICE=D3Q27 -DMPI -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-O1 -g -fcheck=all -cpp -DPRC=4 -DMPI -DVERBOSE -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi	

gnu-real:
	$(MAKE) CC=gcc \
	FC=gfortran \
	CFLAGS="-Ofast -I $(CURDIR) -c" \
	FFLAGS="-Ofast -cpp -DPRC=4 -DLATTICE=D3Q27 -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-Ofast -cpp -DPRC=4 -DVERBOSE -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

gnu-mpi-real:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-Ofast -I $(CURDIR) -c" \
	FFLAGS="-Ofast -cpp -DPRC=4 -DLATTICE=D3Q27 -DMPI -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-Ofast -cpp -DMPI -DPRC=4 -DVERBOSE -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi

gnu-double:
	$(MAKE) CC=gcc \
	FC=gfortran \
	CFLAGS="-Ofast -I $(CURDIR) -c" \
	FFLAGS="-Ofast -cpp -DPRC=8 -DLATTICE=D3Q27 -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-Ofast -cpp -DPRC=8 -DVERBOSE -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

gnu-mpi-double:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-Ofast -I $(CURDIR) -c" \
	FFLAGS="-Ofast -cpp -DPRC=8 -DLATTICE=D3Q27 -DMPI -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-Ofast -cpp -DMPI -DPRC=8 -DVERBOSE -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi

seq:vars_single_3D_module.o \
	boundary_cds_3D_module.o print_module.o recursive_TSLB3D.o
	$(FC) $(LDFLAGS) $(EX) vars_single_3D_module.o \
	boundary_cds_3D_module.o print_module.o recursive_TSLB3D.o
#	mv $(EX) $(EXE)

mpi:vars_single_3D_module.o \
	boundary_cds_3D_module.o print_module.o recursive_TSLB3D.o
	$(FC) $(LDFLAGS) $(EX) vars_single_3D_module.o \
	boundary_cds_3D_module.o print_module.o recursive_TSLB3D.o
#	mv $(EXP) $(EXEP)

get_mem.o:get_mem.c
	$(CC) $(CFLAGS) get_mem.c

get_ram.o:get_ram.c
	$(CC) $(CFLAGS) get_ram.c

vars_single_3D_module.o:vars_single_3D_module.f90
	$(FC) $(FFLAGS) vars_single_3D_module.f90

boundary_cds_3D_module.o:boundary_cds_3D_module.f90
	$(FC) $(FFLAGS) boundary_cds_3D_module.f90

print_module.o:print_module.f90
	$(FC) $(FFLAGS) print_module.f90

recursive_TSLB3D.o:recursive_TSLB3D.f90
	$(FC) $(FFLAGS) recursive_TSLB3D.f90

clean-all:
	rm -rf *.mod *.o *.x 

clean:
	rm -rf *.mod *.o


