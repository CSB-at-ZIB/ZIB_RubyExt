pwd = $(shell pwd)
bpath =$(shell dirname $(pwd))
rpath = $(bpath)/lib

CXX = g++
CXXFLAGS = -I../include -O3
# CXXFLAGS = -g

FC = gfortran
FFLAGS = -O3

OBJS = sbml2adolc.o SbmlModel.o 

all : ../lib/libadolc.so ../lib/libsbml.so sbml2adolc chk_ydot

../lib/libadolc.so :
	cd .. && \
	/bin/rm -rf ADOL-C-2.5.2 && \
	tar xf libADOLC-2.5.2-src.tar.gz \
	cd ADOL-C-2.5.2 && \
	./configure --prefix=$(bpath) --enable-atrig-erf --enable-tserrno --disable-static && \
	$(MAKE) install


../lib/libsbml.so : 
	cd .. && \
	/bin/rm -rf libsbml-4.3.1 && \
	tar xf libSBML-4.3.1-src.tar.gz && \
	cd libsbml-4.3.1 && \
	./configure --prefix=$(bpath) --disable-static && \
	$(MAKE) install

sbml2adolc : $(OBJS)
	$(CXX) -I../include -o sbml2adolc $(OBJS) -L../lib -lsbml -ladolc -Wl,-rpath,$(rpath)
	cp sbml2adolc YDOT_LIMEXcpp_TEMPLATE.txt ../../bin/

SbmlModel.o : SbmlModel.cpp SbmlModel.h myTypes.h

chk_ydot : chk_ydot.o sbml_mathfuns.o ydot_LIMEXcpp.o ydot_LIMEXcpp.h
	$(CXX) -o chk_ydot chk_ydot.o sbml_mathfuns.o ydot_LIMEXcpp.o

clean:
	# /bin/rm -rf ../lib/libadolc.*
	/bin/rm -rf ../ADOL-C-2.5.2
	# /bin/rm -rf ../lib/libsbml.*
	/bin/rm -rf ../libsbml-4.3.1 
	/bin/rm -f *.o sbml2adolc chk_ydot

mclean:
	/bin/rm -rf ../lib/libadolc.*
	/bin/rm -rf ../ADOL-C-2.5.5
	/bin/rm -rf ../lib/libsbml.*
	/bin/rm -rf ../libsbml-4.3.1 
	/bin/rm -f *.o sbml2adolc chk_ydot

