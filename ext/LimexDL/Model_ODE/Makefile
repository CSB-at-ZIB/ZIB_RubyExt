FC = gfortran
FFLAGS = -fexceptions -O3 -fPIC

OBJS = ydot_LIMEX.o sbml_mathfuns.o

all: libODEydot.so

libODEydot.so: $(OBJS) ydot_LIMEX.h
	$(FC) -shared $(OBJS) -o $@ -s 
	# cp $@ ..

clean:
	/bin/rm -f *.so *.o

mclean:
	/bin/rm -f *.so *.o
