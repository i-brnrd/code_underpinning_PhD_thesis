#makefile 1
#hiya : hello.o hello_isla.o
#	#gfortran -o hiya hello.o hello_isla.o

#hello.o : hello.f90
#	#gfortran  hello.f90
#hello_isla.o : hello_isla.f90
#	#gfortran hello_isla.f90
#FFLAGS = -ffpe-trap=invalid,zero -g -fcheck=all -fbacktrace 

objects =search_bisec.o\
	get_dim.o\
	load.o\
	load_spec.o\
	get_cdf.o\
	basic_nogrid.o
#Comp= gfortran

nogrid : $(objects)
	gfortran -o nogrid $(objects)

%.o : %.f90
	gfortran -c -g -fcheck=all -fbacktrace $<
