.SUFFIXES: .f .F .o .a  .f90 .F90
.PHONY: fdf
include ../arch.make

PROGRAM = main

OBJS = $(PROGRAM).o constants.o mpi.o ed_basis.o \
	   numeric_utils.o ed_hamiltonian.o diag.o ed_lattice.o ed_params.o

default: $(PROGRAM)
##################
### FDF Module ###
##################
FDF=fdf/libfdf.a

LIBS+=$(FDF)
FFLAGS+=-I./fdf

####################
### Dependencies ###
####################
$(OBJS): $(FDF)
$(PROGRAM).o: mpi.o ed_params.o ed_basis.o diag.o ed_lattice.o
diag.o: ed_hamiltonian.o
ed_hamiltonian.o: ed_basis.o ed_lattice.o
ed_lattice.o: ed_params.o
ed_basis.o: numeric_utils.o
ed_params.o: numeric_utils.o

$(FDF): 
	cd fdf ; $(MAKE);

$(PROGRAM): $(OBJS)
	$(FC) -o $(PROGRAM).x $(LDFLAGS) $(OBJS) $(LIBS) 

clean:
	rm -f *.x *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)

%.o:%.mod
