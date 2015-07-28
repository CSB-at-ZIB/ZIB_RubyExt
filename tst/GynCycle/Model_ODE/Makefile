FC = gfortran
FINC = -I.
FFLAGS = $(FINC) -fexceptions -O3 -fPIC

CXX = g++
CXXINC = -I.
CXXLIB = -ladolc
CXXFLAGS = $(CXXINC) -fexceptions -O3 -fPIC

OBJS = ydot_LIMEXcpp.o
# OBJS = ydot_LIMEX.o sbml_mathfuns.o

ALL: libODEydot.so

libODEydot.so: $(OBJS) 
	# $(FC) -shared $(OBJS) -o $@ -s
	$(CXX) -shared $(OBJS) -o $@ -s $(CXXLIB)
	cp $@ ..

%.so : %.xml
	./sbml2adolc $< >ydot_LIMEXcpp.cpp
	$(CXX) $(CXXFLAGS) -c ydot_LIMEXcpp.cpp
	$(CXX) -shared ydot_LIMEXcpp.o -o $@ -s $(CXXLIB)
	cp $@ ../libODEydot.so

%.f.so : %.xml
	./sbml2fortran $< >ydot_LIMEX.f
	$(FC) $(FFLAGS) -c ydot_LIMEX.f
	$(FC) $(FFLAGS) -c sbml_mathfuns.f
	$(FC) -shared ydot_LIMEX.o sbml_mathfuns.o -o $@ -s
	cp $@ ../libODEydot.so

clean:
	/bin/rm -f *.so *.o
