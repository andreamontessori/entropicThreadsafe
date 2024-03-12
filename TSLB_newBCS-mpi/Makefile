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
EXP=main.x
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
	@echo "nvfortran-noacc                    "
	@echo "nvfortran-mpi-noacc                "
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
	FFLAGS="-O3 -cpp -acc -gpu=cc75 -Minfo=accel -DPRC=4 -DLATTICE=D3Q27 -I $(CURDIR) -c" \
	LDFLAGS="-O3 -cpp -acc -gpu=cc75 -Minfo=accel -DPRC=4 -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

nvfortran-debug:
	$(MAKE) CC=nvcc \
	FC=nvfortran \
	CFLAGS="-O0 -g -I $(CURDIR) -c" \
	FFLAGS="-O0 -g -Minfo=all -cpp -acc -gpu=cc75 -DPRC=4 -DLATTICE=D3Q27 -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-O0 -g -cpp -acc -gpu=cc75 -DPRC=4 -DVERBOSE -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

nvfortran-mpi:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O3 -I $(CURDIR) -c" \
	FFLAGS="-O3 -cpp -acc -gpu=cc75 -Minfo=accel -DPRC=4 -DLATTICE=D3Q27 -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-O3 -cpp -acc -gpu=cc75 -Minfo=accel -DPRC=4 -DMPI -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi

nvfortran-noacc:
	$(MAKE) CC=nvcc \
	FC=nvfortran \
	CFLAGS="-O3 -I $(CURDIR) -c" \
	FFLAGS="-O3 -cpp -DPRC=4 -DLATTICE=D3Q27 -I $(CURDIR) -c" \
	LDFLAGS="-O3 -cpp -DPRC=4 -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

nvfortran-mpi-noacc:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O1 -I $(CURDIR) -c" \
	FFLAGS="-O1 -cpp -DPRC=4 -DLATTICE=D3Q27 -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-O1 -cpp -DPRC=4 -DMPI -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi

gnu-debug:
	$(MAKE) CC=gcc \
	FC=gfortran \
	CFLAGS="-O1 -g -I $(CURDIR) -c" \
	FFLAGS="-O1 -g -fcheck=all -ffree-line-length-none -cpp -DPRC=4 -DLATTICE=D3Q27 -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-O1 -g -fcheck=all -ffree-line-length-none -cpp -DPRC=4 -DVERBOSE -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

gnu-mpi-debug:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O1 -g -I $(CURDIR) -c" \
	FFLAGS="-O1 -g -fcheck=all -ffree-line-length-none -cpp -DPRC=4 -DLATTICE=D3Q27 -DMPI -DVERBOSE -I $(CURDIR) -c" \
	LDFLAGS="-O1 -g -fcheck=all -ffree-line-length-none -cpp -DPRC=4 -DMPI -DVERBOSE -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi	

gnu-real:
	$(MAKE) CC=gcc \
	FC=gfortran \
	CFLAGS="-Ofast -I $(CURDIR) -c" \
	FFLAGS="-Ofast -ffree-line-length-none -cpp -DPRC=4 -DLATTICE=D3Q27 -I $(CURDIR) -c" \
	LDFLAGS="-Ofast -ffree-line-length-none -cpp -DPRC=4 -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

gnu-mpi-real:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-Ofast -I $(CURDIR) -c" \
	FFLAGS="-Ofast -ffree-line-length-none -cpp -DPRC=4 -DLATTICE=D3Q27 -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-Ofast -ffree-line-length-none -cpp -DMPI -DPRC=4 -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi

gnu-double:
	$(MAKE) CC=gcc \
	FC=gfortran \
	CFLAGS="-Ofast -I $(CURDIR) -c" \
	FFLAGS="-Ofast -ffree-line-length-none -cpp -DPRC=8 -DLATTICE=D3Q27 -I $(CURDIR) -c" \
	LDFLAGS="-Ofast -ffree-line-length-none -cpp -DPRC=8 -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq

gnu-mpi-double:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-Ofast -I $(CURDIR) -c" \
	FFLAGS="-Ofast -ffree-line-length-none -cpp -DPRC=8 -DLATTICE=D3Q27 -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-Ofast -ffree-line-length-none -cpp -DMPI -DPRC=8 -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi

seq:get_mem.o get_ram.o vars_single_3D_module.o \
	mpi_module.o profiling_m.o boundary_cds_3D_module.o lb_kernels_module.o print_module.o recursive_TSbcs_TSLB3.o
	$(FC) $(LDFLAGS) $(EX) get_mem.o get_ram.o vars_single_3D_module.o \
	mpi_module.o profiling_m.o boundary_cds_3D_module.o lb_kernels_module.o print_module.o recursive_TSbcs_TSLB3.o
#	mv $(EX) $(EXE)

mpi:get_mem.o get_ram.o vars_single_3D_module.o \
	mpi_module.o profiling_m.o boundary_cds_3D_module.o lb_kernels_module.o print_module.o recursive_TSbcs_TSLB3.o
	$(FC) $(LDFLAGS) $(EX) get_mem.o get_ram.o vars_single_3D_module.o \
	mpi_module.o profiling_m.o boundary_cds_3D_module.o lb_kernels_module.o print_module.o recursive_TSbcs_TSLB3.o
#	mv $(EXP) $(EXEP)

get_mem.o:get_mem.c
	$(CC) $(CFLAGS) get_mem.c

get_ram.o:get_ram.c
	$(CC) $(CFLAGS) get_ram.c

vars_single_3D_module.o:vars_single_3D_module.f90
	$(FC) $(FFLAGS) vars_single_3D_module.f90
	
mpi_module.o:mpi_module.f90
	$(FC) $(FFLAGS) mpi_module.f90

profiling_m.o: profiling_m.f90
	$(FC) $(FFLAGS) profiling_m.f90

boundary_cds_3D_module.o:boundary_cds_3D_module.f90
	$(FC) $(FFLAGS) boundary_cds_3D_module.f90

lb_kernels_module.o:lb_kernels_module.f90
	$(FC) $(FFLAGS) lb_kernels_module.f90

print_module.o:print_module.f90
	$(FC) $(FFLAGS) print_module.f90

recursive_TSLB3D.o:recursive_TSLB3D.f90
	$(FC) $(FFLAGS) recursive_TSLB3D.f90

recursive_TSbcs_TSLB3.o:recursive_TSbcs_TSLB3.f90
	$(FC) $(FFLAGS) recursive_TSbcs_TSLB3.f90

clean-all:
	rm -rf *.mod *.o *.x *.dat

clean:
	rm -rf *.mod *.o


