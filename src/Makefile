MAKEFLAGS = -s
FC = ifort
#FCFLAGS = -static -static-intel -static-libgcc -mkl -traceback -w -I. -g -fpp -Iutils/
FCFLAGS = -mkl -traceback -w -I. -g -fpp -Iutils/
LDFLAGS = -Lutils/ -lutils
TLIB = 

all: utils/libutils.a InputParser.o MixedState.o MixedStateTransition.o PureState.o PureStateTransition.o  $(TLIB) HenX2Spec 

utils/libutils.a:
	make -C utils

HenX2Spec: utils/libutils.a HenX2Spec.f90 PureStateTransition.o PureState.o MixedState.o MixedStateTransition.o InputParser.o purestatetransition_.mod purestate_.mod mixedstate_.mod mixedstatetransition_.mod inputparser_.mod
	echo 'Building HenX2Spec.f90'
	$(FC) $(FCFLAGS) HenX2Spec.f90 PureStateTransition.o PureState.o MixedState.o MixedStateTransition.o InputParser.o -o HenX2Spec $(LDFLAGS)

InputParser.o inputparser_.mod: InputParser.f90 
	echo 'Building InputParser.f90'
	$(FC) $(FCFLAGS) -c InputParser.f90

MixedState.o mixedstate_.mod: MixedState.f90 purestate_.mod
	echo 'Building MixedState.f90'
	$(FC) $(FCFLAGS) -c MixedState.f90

MixedStateTransition.o mixedstatetransition_.mod: MixedStateTransition.f90 purestate_.mod mixedstate_.mod
	echo 'Building MixedStateTransition.f90'
	$(FC) $(FCFLAGS) -c MixedStateTransition.f90

PureState.o purestate_.mod: PureState.f90 
	echo 'Building PureState.f90'
	$(FC) $(FCFLAGS) -c PureState.f90

PureStateTransition.o purestatetransition_.mod: PureStateTransition.f90 purestate_.mod
	echo 'Building PureStateTransition.f90'
	$(FC) $(FCFLAGS) -c PureStateTransition.f90

clean:
	make -C utils clean
	rm -f *.o *.mod *~ HenX2Spec  2> /dev/null

