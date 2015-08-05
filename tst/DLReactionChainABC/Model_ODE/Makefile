pwd = $(shell pwd)
bpath =$(shell dirname $(pwd))
rpath = $(bpath)/../../pkg/lib

FC = gfortran
FINC = -I.
FFLAGS = $(FINC) -fexceptions -O3 -fPIC

CXX = g++
CXXINC = -I. 
CXXLIB = -ladolc 
CXXINC = -I. -I../../../pkg/include
CXXLIB = -Wl,-init,set_adolc_num_dir -L../../../pkg/lib -ladolc -Wl,-rpath,$(rpath)
CXXFLAGS = $(CXXINC) -fexceptions -O3 -fPIC

OBJS = ydot_LIMEXcpp.o
#OBJS = ydot_LIMEX.o sbml_mathfuns.o

ALL: libODEydot.so

libODEydot.so: $(OBJS) 
	$(CXX) -s -shared $(OBJS) -o $@ $(CXXLIB)
	cp $@ ..

%.so : %.xml
	./sbml2adolc $< >ydot_LIMEXcpp.cpp
	$(CXX) $(CXXFLAGS) -c ydot_LIMEXcpp.cpp
	$(CXX) -s -shared ydot_LIMEXcpp.o -o $@ $(CXXLIB)
	cp $@ ../libODEydot.so

%.f.so : %.xml
	./sbml2fortran $< >ydot_LIMEX.f
	$(FC) $(FFLAGS) -c ydot_LIMEX.f
	$(FC) $(FFLAGS) -c sbml_mathfuns.f
	$(FC) -s -shared ydot_LIMEX.o sbml_mathfuns.o -o $@
	cp $@ ../libODEydot.so

clean:
	/bin/rm -f *.so *.o
