FORTRAN = gfortran

all: data
	${FORTRAN} -g *.f -o GM_cryocooler.o
data:
	${FORTRAN} -c data.f
clean:
	rm -f *.o *.mod
	rm -f Pressure.txt
	rm -f Flux.txt
	rm -f GMout.txt
	rm -f *.png
	rm -f *~
