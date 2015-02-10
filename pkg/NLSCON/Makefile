FC = gfortran
FFLAGS = -fexceptions -O3 -fPIC
OBJS  =  nlscon.o wnorm.o linalg_nlscon.o zibconst.o zibsec.o zibmon.o

all : libNLSCON.so

libNLSCON.so: $(OBJS)
	$(FC) -shared $(OBJS) -o $@ -s
	# cp $@ ..

nlscon: main_nlscon.o $(OBJS)
	$(FC) $(DFLAGS) main_nlscon.o $(OBJS) -o nlscon

clean:
	/bin/rm -f *.o *.so nlscon

mclean:
	/bin/rm -f *.o *.so nlscon
