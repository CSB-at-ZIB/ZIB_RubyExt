#
#  Makefile for LIMEX version 4.3A
#  -------------------------------
#
#  FORTRAN sources:
#
#     LIMEX4_3A.f, LIMEX4_3_Dense.f,
#
#  together with one of the sample programs:
#
#     LIMEX4_3A_Example_Small.f or LIMEX4_3A_Example_Big.f
# 
#  furthermore, if BLAS and LAPACK libraries are not available:
#
#     LIMEX4_3_Auxiliaries.f
#  
LOBJ  =  LIMDHERM.o \
         LIMEX4_3A.o LIMEX4_3_Dense.o \
         LIMEX4_3_Auxiliaries.o
#
#  FORTRAN compiler options:
#
FC = gfortran
FFLAGS = -fexceptions -O3 -fPIC
#
all : libLIMEX4_3A.so
#
libLIMEX4_3A.so : $(LOBJ)
	$(FC) -shared $(LOBJ) -o $@ -s 
	# cp $@ ..
#
#  Linked without BLAS and LAPACK libraries:
#
limex : $(LOBJ) LIMEX4_3A_Example_Small.o
	$(FC) -o $@ $(LOBJ) LIMEX4_3A_Example_Small.o
#
#  Linked with BLAS and LAPACK libraries (e.g. myblas.a and mylapack.a):
#
#limex:  $(LOBJ)
#	f77 -o limex $(LOBJ) myblas.a mylapack.a
#
#  LIMEX4_A.f depends from LIMEX4_3_Size_Definitions.h:
#
LIMEX4_3A.o : LIMEX4_3A.f   LIMEX4_3_Size_Definitions.h
LIMDHERM.o : LIMDHERM.f   LIMEX4_3_Size_Definitions.h

clean :
	/bin/rm -f *.o *.so

mclean :
	/bin/rm -f *.o *.so
