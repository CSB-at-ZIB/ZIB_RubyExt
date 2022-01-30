FC = gfortran
FFLAGS = -fexceptions -O3 -fPIC
OBJS  =  $(patsubst %.f,%.o,$(shell ls *.f))

all : libDIERCKX.so

libDIERCKX.so: $(OBJS)
	$(FC) -shared $(OBJS) -o $@ -s
	# cp $@ ..

#nlscon: main_nlscon.o $(OBJS)
#	$(FC) $(FFLAGS) main_nlscon.o $(OBJS) -o nlscon

clean:
	/bin/rm -f *.o *.so 

mclean:
	/bin/rm -f *.o *.so 
