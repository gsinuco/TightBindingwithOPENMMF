GF=gfortran

###################################
# MAKE CPP EXAMPLES
###################################

all: zvode

zvode: zvode_call.f
	gfortran -c zgbfa.f
	gfortran -c zgbsl.f
	gfortran -c zgefa.f
	gfortran -c zgesl.f
	gfortran -c zvode.f
	gfortran -o zvode zgbfa.o zgbsl.o zgefa.o zgesl.o zvode.o zvode_call.f -lblas


############################
# CLEAN
############################

clean:
	rm *.o zvode

