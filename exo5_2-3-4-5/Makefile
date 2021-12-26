##########################################
# Makefile                               #
# Makefile for the code developed in TP1 #
#                                        #
# T. Dufaud                              #
##########################################
################################
# Variables for this makefile
################################
# 
CC=gcc

# 
# -- Compiler Option
#
OPTC=-O3 -fomit-frame-pointer -fPIC -mavx -DAdd_ -DF77_INTEGER=int -DStringSunStyle
# These macros identifies how FORTRAN routines wil be called : 
#  Add_     : name of C functions to be in all lower case and to have an underscore postfixed it (Suns, Intel compilers expect this).
#  NoChange : the Fortran compiler expects the name of C functions to be in all lower case (IBM RS6K compilers do this). 
#  UpCase   : the Fortran compiler expects the name of C functions to be in all upcase (Cray compilers expect this). 
#  Add__    : the Fortran compiler in use is f2c, a Fortran to C converter.
#  F77_INTEGER=int : 
#  StringSunStyle :
#
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src

#
# -- librairies
LIBS=-llapacke -lcblas -lm

# -- Include directories
INCLBLASLAPACK= -I /usr/include

INCL= -I $(TPDIR)/include $(INCLBLASLAPACK) 
#
#################################################################
# makefile
############
#
OBJENV= tp_env.o
OBJTP2ITER= lib_poisson1D.o tp2_poisson1D_iter.o
OBJTP2DIRECT= lib_poisson1D.o tp2_poisson1D_direct.o
OBJTP2DGBMV= lib_poisson1D.o tp2_poisson1D_dgbmv.o
OBJTP2LUTRID= lib_poisson1D.o tp2_poisson1D_lutrid.o
#

all: bin/tp_testenv bin/tp2poisson1D_iter bin/tp2poisson1D_direct bin/tp2poisson1D_dgbmv bin/tp2poisson1D_lutrid

testenv: bin/tp_testenv

tp2poisson1D_iter: bin/tp2poisson1D_iter

tp2poisson1D_direct: bin/tp2poisson1D_direct

tp2poisson1D_dgbmv: bin/tp2poisson1D_dgbmv

tp2poisson1D_lutrid: bin/tp2poisson1D_lutrid

tp_env.o: $(TPDIRSRC)/tp_env.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_env.c 

lib_poisson1D.o: $(TPDIRSRC)/lib_poisson1D.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/lib_poisson1D.c 

tp2_poisson1D_iter.o: $(TPDIRSRC)/tp2_poisson1D_iter.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_iter.c  

tp2_poisson1D_direct.o: $(TPDIRSRC)/tp2_poisson1D_direct.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_direct.c  

tp2_poisson1D_dgbmv.o: $(TPDIRSRC)/tp2_poisson1D_dgbmv.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_dgbmv.c

tp2_poisson1D_lutrid.o: $(TPDIRSRC)/tp2_poisson1D_lutrid.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_lutrid.c  

bin/tp_testenv: $(OBJENV) 
	$(CC) -o bin/tp_testenv $(OPTC) $(OBJENV) $(LIBS)

bin/tp2poisson1D_iter: $(OBJTP2ITER)
	$(CC) -o bin/tp2poisson1D_iter $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tp2poisson1D_direct: $(OBJTP2DIRECT)
	$(CC) -o bin/tp2poisson1D_direct $(OPTC) $(OBJTP2DIRECT) $(LIBS)

bin/tp2poisson1D_dgbmv: $(OBJTP2DGBMV)
	$(CC) -o bin/tp2poisson1D_dgbmv $(OPTC) $(OBJTP2DGBMV) $(LIBS)

bin/tp2poisson1D_lutrid: $(OBJTP2LUTRID)
	$(CC) -o bin/tp2poisson1D_lutrid $(OPTC) $(OBJTP2LUTRID) $(LIBS)

run_testenv:
	bin/tp_testenv

run_tp2poisson1D_iter:
	bin/tp2poisson1D_iter

run_tp2poisson1D_direct:
	bin/tp2poisson1D_direct

run_tp2poisson1D_dgbmv:
	bin/tp2poisson1D_dgbmv

run_tp2poisson1D_lutrid:
	bin/tp2poisson1D_lutrid

clean:
	rm *.o *.dat bin/*
