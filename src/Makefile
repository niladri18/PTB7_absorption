version = 05-21-2014


FC = ifort
FCFLAGS = -O1 -ipo
LIB = -L/usr/local/lib -mkl




PROGRAM = sci

all = $(PROGRAM)

all: param.o read.o make_huckel.o hartree.o sci.o
	$(FC) $(FCFLAGS) -o sci param.o read.o make_huckel.o hartree.o sci.o ${LIB}

param.o: param.f90
	$(FC) $(FCFLAGS) -c param.f90 $(LIB)
read.o: read.f90
	$(FC) $(FCFLAGS) -c read.f90 $(LIB)
make_huckel.o: make_huckel.f90
	$(FC) $(FCFLAGS) -c make_huckel.f90 $(LIB)
hartree.o: hartree.f90
	$(FC) $(FCFLAGS) -c hartree.f90 $(LIB)
sci.o: sci.f90
	 $(FC) $(FCFLAGS) -c sci.f90  $(LIB)
#pirg-rtc.o: pirg-rtc.f90
#
##$(FC) $(FCFLAGS) -c pirg-rtc.f90 $(LIB)






clean: 
	rm -f pirg setup.o variable.o inv.o lapack.o utils.o energy.o utype.o pirg10.o *.o, *.mod,*~

