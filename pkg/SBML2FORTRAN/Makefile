pwd = $(shell pwd)
bpath =$(shell dirname $(pwd))
rpath = $(bpath)/lib

CXX = g++
CXXFLAGS = -I../include -O3
# CXXFLAGS = -g

FC = gfortran
FFLAGS = -O3

OBJS = sbml2fortran.o SbmlModel.o SbmlJacobian.o \
       Parser.o Expression.o ExprNode.o

all : ../lib/libsbml.so sbml2fortran chk_ydot

../lib/libsbml.so : 
	cd .. && \
	/bin/rm -rf libsbml-4.3.1 && \
	tar xf libSBML-4.3.1-src.tar.gz && \
	cd libsbml-4.3.1 && \
	./configure --prefix=$(bpath) --disable-static && \
	$(MAKE) install

sbml2fortran : $(OBJS)
	$(CXX) -I../include -o sbml2fortran $(OBJS) -L../lib -lsbml -Wl,-rpath $(rpath)
	cp sbml2fortran YDOT_LIMEX_TEMPLATE.txt ../../bin/

SbmlModel.o : SbmlModel.cpp SbmlModel.h myTypes.h

SbmlJacobian.o : SbmlJacobian.cpp SbmlJacobian.h myTypes.h

Parser.o : Parser.cpp Parser.h Expression.h myTypes.h

Expression.o : Expression.cpp Expression.h myTypes.h

ExprNode.o : ExprNode.cpp ExprNode.h Expression.h myTypes.h

chk_ydot : chk_ydot.o sbml_mathfuns.o ydot_LIMEX.o ydot_LIMEX.h
	$(CXX) -o chk_ydot chk_ydot.o sbml_mathfuns.o ydot_LIMEX.o

clean:
	# /bin/rm -rf ../lib/libsbml.*
	/bin/rm -rf ../libsbml-4.3.1 
	/bin/rm -f *.o sbml2fortran chk_ydot

mclean:
	/bin/rm -rf ../lib/libsbml.*
	/bin/rm -rf ../libsbml-4.3.1 
	/bin/rm -f *.o sbml2fortran chk_ydot
