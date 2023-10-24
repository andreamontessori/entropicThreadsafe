BINROOT := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))

FC=nvfortran
CUDAFLAGS = -cuda -gpu=cc70 -O0 -g -Mbounds -Mchkptr -Mchkstk
CUDAFLAGS = -cuda -gpu=cc70,keepptx -O0 -g
CUDAFLAGS = -cuda -fast -gpu=cc80 -DMYDIMESION=128 -DTILE1=128 -DTILE2=1 -DTILE3=1
CUDAFLAGS = -gpu=cc75,lineinfo,ptxinfo -cpp -cuda -O3
CUDAFLAGSREG = -gpu=cc75,cuda11.0,lineinfo,ptxinfo,maxregcount:128 -cpp -cuda -O3
LDFLAGS = -gpu=cc75,cuda11.0,lineinfo,ptxinfo,maxregcount:128 -cpp -cuda -O3 -o 

def:	all

all:
	@echo "Error - you should specify a target machine!"
	@echo "Possible choices are:              "
	@echo "                                   "
	@echo "2D                                 "
	@echo "3D                                 "
	@echo "                                   "
	@echo "Please examine Makefile for further details "

help: all
	
2D:
	$(MAKE) FC=nvfortran \
	EX=main2D.x \
	TYPE="seq" \
	BINROOT=$(BINROOT) seq2d
	
3D:
	$(MAKE) FC=nvfortran \
	TYPE=mpi \
	EX=main3D.x \
	BINROOT=$(BINROOT) seq3d

seq2d:vars_single_2D.o print_mod.o recTSLB_thirdorder.o 
	$(FC) $(LDFLAGS) $(EX) vars_single_2D.o print_mod.o \
	recTSLB_thirdorder.o

seq3d:vars_single_3D.o print_mod.o recTSLB_thirdorder3D.o
	$(FC) $(LDFLAGS) $(EX) vars_single_3D.o print_mod.o \
	recTSLB_thirdorder3D.o

vars_single_2D.o: Makefile vars_single_2D.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c vars_single_2D.f90

vars_single_3D.o: Makefile vars_single_3D.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c vars_single_3D.f90

print_mod.o: Makefile print_mod.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c print_mod.f90

recTSLB_thirdorder.o: Makefile recTSLB_thirdorder.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c recTSLB_thirdorder.f90

recTSLB_thirdorder3D.o: Makefile recTSLB_thirdorder3D.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c recTSLB_thirdorder3D.f90

clean:
	@rm -rf *.x *.o *.mod *.ptx
