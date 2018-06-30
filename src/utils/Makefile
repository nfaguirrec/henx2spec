MAKEFLAGS = -s
FC = ifort
FCFLAGS = -traceback -w -I. -g -fpp
LDFLAGS = 
TLIB = libutils.a

all: Grid.o Math.o Matrix.o Morse.o NFunction.o NIntegration.o ThrularNumerovMethod.o UnitsConverter.o  $(TLIB) 

Grid.o grid_.mod: Grid.f90 
	echo 'Building Grid.f90'
	$(FC) $(FCFLAGS) -c Grid.f90

Math.o math_.mod: Math.f90 
	echo 'Building Math.f90'
	$(FC) $(FCFLAGS) -c Math.f90

Matrix.o matrix_.mod: Matrix.f90 
	echo 'Building Matrix.f90'
	$(FC) $(FCFLAGS) -c Matrix.f90

Morse.o morse_.mod: Morse.f90 unitsconverter_.mod nfunction_.mod thrularnumerovmethod_.mod grid_.mod
	echo 'Building Morse.f90'
	$(FC) $(FCFLAGS) -c Morse.f90

NFunction.o nfunction_.mod: NFunction.f90 grid_.mod
	echo 'Building NFunction.f90'
	$(FC) $(FCFLAGS) -c NFunction.f90

NIntegration.o nintegration_.mod: NIntegration.f90 nfunction_.mod grid_.mod
	echo 'Building NIntegration.f90'
	$(FC) $(FCFLAGS) -c NIntegration.f90

ThrularNumerovMethod.o thrularnumerovmethod_.mod: ThrularNumerovMethod.f90 nfunction_.mod grid_.mod
	echo 'Building ThrularNumerovMethod.f90'
	$(FC) $(FCFLAGS) -c ThrularNumerovMethod.f90

UnitsConverter.o unitsconverter_.mod: UnitsConverter.f90 
	echo 'Building UnitsConverter.f90'
	$(FC) $(FCFLAGS) -c UnitsConverter.f90

$(TLIB): Grid.o Math.o Matrix.o Morse.o NFunction.o NIntegration.o ThrularNumerovMethod.o UnitsConverter.o 
	echo 'Building library $(TLIB)'
	ar rvs $(TLIB) Grid.o Math.o Matrix.o Morse.o NFunction.o NIntegration.o ThrularNumerovMethod.o UnitsConverter.o  > /dev/null 2> /dev/null

clean:
	rm -f *.o *.mod *~ $(TLIB) 2> /dev/null

