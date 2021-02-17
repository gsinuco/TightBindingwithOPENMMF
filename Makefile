#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:../../lib/"
TOPSRCDIR = ./
include $(TOPSRCDIR)/make.inc

#SET MULTIMODEFLOQUET LIBRARY PATH
MMFLIB = ./lib/
#SET MULTIMODEFLOQUET INCLUDE PATH
MMFINC = ./include/

###################################
# MAKE CPP EXAMPLES
###################################

ifndef BUILD_MKL
BUILD_MKL_ = no
all : ManybodyHubbard 
endif

ifdef BUILD_MKL
BUILD_MKL_ = yes
all: ManybodyHubbard
endif



#SET MULTIMODEFLOQUET LIBRARY PATH
MMFLIB = ./lib/
#SET MULTIMODEFLOQUET INCLUDE PATH
MMFINC = ./include/

###################################
# MAKE CPP EXAMPLES
###################################

all: ManybodyHubbard ManybodyHubbard_RHS #

ManybodyHubbard: subset.f90 HilbertDimension.f90 main_MBH_Bosons.f90
	$(GF) -c bvec.f90
	$(GF) -c subset.f90
	$(GF) -c HilbertDimension.f90
	$(GF) -c HubbardHamiltonian.f90
	$(GF) -c fftw_1D.f90  -I/usr/include/ -O3 -llapack -lblas -g -lstdc++ -lfftw3	
	$(GF) -c tight_binding_boson_Hamiltonian.f90 -I./include
	$(GF) -c Runge-Kutta4thorder_SIMPLE.f90
	$(GF) -c LapackEigenValues.f90
	$(GF) -o ManybodyHubbard_Bosons  LapackEigenValues.o tight_binding_boson_Hamiltonian.o Runge-Kutta4thorder_SIMPLE.o fftw_1D.o bvec.o subset.o HilbertDimension.o HubbardHamiltonian.o main_MBH_Bosons.f90 -I/usr/include -I$(MMFINC) -L$(MMFLIB) $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(MKLFLAGS)



ManybodyHubbard_RHS: subset.f90 HilbertDimension.f90 main_MBH_Bosons.f90
	$(GF) -c bvec.f90
	$(GF) -c subset.f90
	$(GF) -c HilbertDimension.f90
	$(GF) -c HubbardHamiltonian.f90
	$(GF) -c fftw_1D.f90  -I/usr/include/ -O3 -llapack -lblas -g -lstdc++ -lfftw3	
	$(GF) -c tight_binding_boson_Hamiltonian.f90 -I./include
	$(GF) -c Runge-Kutta4thorder_RHS.f90
	$(GF) -c LapackEigenValues.f90
	$(GF) -o ManybodyHubbard_Bosons  LapackEigenValues.o tight_binding_boson_Hamiltonian.o Runge-Kutta4thorder_RHS.o fftw_1D.o bvec.o subset.o HilbertDimension.o HubbardHamiltonian.o main_MBH_Bosons.f90 -I/usr/include -I$(MMFINC) -L$(MMFLIB) $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(MKLFLAGS)


############################
# CLEAN
############################

clean:
	rm *.o ManybodyHubbard_Fermions ManybodyHubbard_Bosons

