FC = gfortran
FFLAGS = -fexceptions -O3 -fPIC
FINC = -I.

CXX = g++
CXXFLAGS = -fexceptions -O3 -fPIC
CXXINC = -I.

OBJS = ydot_LIMEXcpp.o
# OBJS = ydot_LIMEX.o sbml_mathfuns.o

ALL: libODEydot.so

libODEydot.so: $(OBJS) 
	# $(FC) -shared $(OBJS) -o $@ -s
	$(CXX) -shared $(OBJS) -o $@ -s
	cp $@ ..

%.so : %.xml
	./sbml2adolc $< >ydot_LIMEXcpp.cpp
	$(CXX) $(CXXFLAGS) $(CXXINC) -c ydot_LIMEXcpp.cpp
	$(CXX) -shared ydot_LIMEXcpp.o -o $@ -s
	cp $@ ../libODEydot.so

# %.so : %.xml
# 	./sbml2fortran $< >ydot_LIMEX.f
# 	$(FC) $(FFLAGS) $(FINC) -c ydot_LIMEX.f
# 	$(FC) $(FFLAGS) $(FINC) -c sbml_mathfuns.f
# 	$(FC) -shared ydot_LIMEX.o sbml_mathfuns.o -o $@ -s
# 	cp $@ ../libODEydot.so

clean:
	/bin/rm -f *.so *.o
