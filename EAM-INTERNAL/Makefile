PROG =	bh_v07

SRCS =	parameters.f90 module_bh_mod.f90 module_function_par.f90  module_histo.f90\
        main_bh.f90 read_bh.f90 initialization.f90\
        drive_bh.f90 gauss.f90 ranmar_bh.f  initial_configuration.f90 \
	minimization_rgl.f90 l-bfgs-b.f90 \
	comparison.f90 bigvoi.f90 force_rgl.f90\
	struttura.f90 cna.f90 permutazioni.f90 calcola_histo.f90 \
	output_bh_glo.f90 output_bh_lista.f90 \
	move_ball.f90 move_single.f90 move_shake.f90 move_bonds.f90 \
	move_shell.f90 move_exchange.f90 move_highenergyatoms.f90 move_brownian.f90 move_brownian_surf.f90\
	move_mirror.f90\
	minima_collecting.f90 cont_neigh.f90 standard_bh.f90\
	force_pres.f90 autov.f90\
	grouping.f90\


OBJS =	parameters.o module_bh_mod.o module_function_par.o module_histo.o\
        main_bh.o read_bh.o initialization.o\
        drive_bh.o gauss.o ranmar_bh.o initial_configuration.o \
	minimization_rgl.o l-bfgs-b.o \
	comparison.o bigvoi.o force_rgl.o\
	struttura.o cna.o permutazioni.o calcola_histo.o \
	output_bh_glo.o output_bh_lista.o \
	move_ball.o move_single.o move_shake.o move_bonds.o \
	move_shell.o  move_exchange.o move_highenergyatoms.o move_brownian.o move_brownian_surf.o\
	move_mirror.o\
	minima_collecting.o cont_neigh.o standard_bh.o\
	force_pres.o autov.o\
	grouping.o\

#        buffer.o

#LIBS = 

#FC=pgf90
#FC=gfortran
#FFLAGS=-static-libcxa -Xlinker -Bstatic -lsvml\
#-L/opt/intel/mkl/lib/32/ -lmkl_lapack -lmkl_p4 \
#-lguide -Vaxlib -Xlinker -Bdynamic -lpthread
#FC=ifort
#FC= pgf90
FC=gfortran -O3
FFLAGS= -std=legacy
F90=gfortran -O3
#F90=pgf90
#F90 =pgf90
#F90= ifort  -O1
F90FLAGS= -std=legacy
#LDFLAGS = -Wl,-R/usr/lib/gcc-lib/i486-linux/3.3.5 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c  $<

