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
	FFLAGS="-acc -Mpreprocess -Minfo=all -ta=tesla -O3 -I $(CURDIR) -c" \
	LDFLAGS="-acc -Mpreprocess -Minfo=all -ta=tesla -O3 -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

nvfortran-mpi:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O3 -I $(CURDIR) -c" \
	FFLAGS="-acc -Mpreprocess -Minfo=all -ta=tesla -O3 -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-acc -Mpreprocess -Minfo=all -ta=tesla -O3 -DMPI -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi

seq:get_mem.o vars_single_3D_module.o mpi_module.o \
	print_module.o boundary_cds_3D_module.o recursive_TSLB3D.o
	$(FC) $(LDFLAGS) $(EX) get_mem.o vars_single_3D_module.o \
	mpi_module.o print_module.o boundary_cds_3D_module.o recursive_TSLB3D.o
#	mv $(EX) $(EXE)

mpi:get_mem.o vars_single_3D_module.o mpi_module.o \
	print_module.o boundary_cds_3D_module.o recursive_TSLB3D.o
	$(FC) $(LDFLAGS) $(EX) get_mem.o vars_single_3D_module.o \
	mpi_module.o print_module.o boundary_cds_3D_module.o recursive_TSLB3D.o
#	mv $(EXP) $(EXEP)

get_mem.o:get_mem.c
	$(CC) $(CFLAGS) get_mem.c

vars_single_3D_module.o:vars_single_3D_module.f90
	$(FC) $(FFLAGS) vars_single_3D_module.f90

mpi_module.o:mpi_module.f90
	$(FC) $(FFLAGS) mpi_module.f90

print_module.o:print_module.f90
	$(FC) $(FFLAGS) print_module.f90

boundary_cds_3D_module.o:boundary_cds_3D_module.f90
	$(FC) $(FFLAGS) boundary_cds_3D_module.f90

recursive_TSLB3D.o:recursive_TSLB3D.f90
	$(FC) $(FFLAGS) recursive_TSLB3D.f90

clean:
	rm -rf *.mod *.o main*

