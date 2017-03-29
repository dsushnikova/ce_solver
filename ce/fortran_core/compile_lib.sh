#gfortran  -g -pg  -fprofile-arcs -ftest-coverage putstrmodule.F90 dispmodule.f90 fast_direct.f90 test.f90 -llapack -lblas -o main
#ifort -w -O0  -i4 putstrmodule.F90 dispmodule.f90 fast_direct.f90 test.f90 -mkl -g -pg -p -o main
gfortran  -c -O3 sparsekit.f90 dlauc1.f dgeqpw.f dgeqpc.f dtrrnk.f dlasmx.f dtrqxc.f dgeqpx.f dtrqpx.f dgeqpb.f dsort.f blassm.f matvec.f unary.f formats.f putstrmodule.F90 dispmodule.f90
#ifort  -c -O3 dlauc1.f dgeqpw.f dgeqpc.f dtrrnk.f dlasmx.f dtrqxc.f dgeqpx.f dtrqpx.f dgeqpb.f dsort.f blassm.f unary.f formats.f putstrmodule.F90 dispmodule.f90 
ar rc  mylib.a dlauc1.o dgeqpw.o dgeqpc.o dtrrnk.o dlasmx.o dtrqxc.o dgeqpx.o dtrqpx.o dgeqpb.o     dsort.o blassm.o unary.o formats.o putstrmodule.o dispmodule.o
