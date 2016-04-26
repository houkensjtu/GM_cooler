all: data
	gfortran -free -g *.f -o GM_cryocooler.o
data:
	gfortran -free -c data.f
clean:
	rm -f *.o *.mod
	rm -f GM.txt
	rm -f GMout.txt
