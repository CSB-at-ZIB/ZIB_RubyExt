pwd = $(shell pwd)
bpath = $(shell dirname $(pwd))/../../pkg
pcpath = $(bpath)/lib64/pkgconfig:$(bpath)/lib/pkgconfig
pkg_config_add = $(shell PKG_CONFIG_PATH=$(pcpath) pkg-config --libs adolc)

FC = gfortran
FINC = -I.
FFLAGS = $(FINC) -fexceptions -O3 -fPIC

CXX = g++
CXXINC = -I. -I../../../pkg/include
CXXFLAGS = $(CXXINC) -fexceptions -O3 -fPIC
CXXLIB = -Wl,-init,set_adolc_num_dir $(pkg_config_add)

OBJS = ydot_LIMEXcpp.o
# OBJS = ydot_LIMEX.o sbml_mathfuns.o

ALL: libODEydot.so

libODEydot.so: $(OBJS)
	# $(FC) -shared $(OBJS) -o $@ -s 
	$(CXX) -shared -s $(OBJS) -o $@ $(CXXLIB)
	cp $@ ..

%.so : %.xml
	./sbml2adolc $< >ydot_LIMEXcpp.cpp
	$(CXX) $(CXXFLAGS) -c ydot_LIMEXcpp.cpp
	$(CXX) -shared -s ydot_LIMEXcpp.o -o $@ $(CXXLIB)
	cp $@ ../libODEydot.so

%.f.so : %.xml
	./sbml2fortran $< >ydot_LIMEX.f
	$(FC) $(FFLAGS) -c ydot_LIMEX.f
	$(FC) $(FFLAGS) -c sbml_mathfuns.f
	$(FC) -shared -s ydot_LIMEX.o sbml_mathfuns.o -o $@
	cp $@ ../libODEydot.so

clean:
	/bin/rm -f *.so *.o
