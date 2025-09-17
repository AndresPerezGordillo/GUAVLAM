#
# Makefile for program 'ShroudedPropeller' 
# Author: Juan D. Colmenares
#

SRCS_f901 = uvlm.f90 \
	   surfaceModule.f90 \
	   wakeModule.f90 \
	   solvers.f90 \
	   tecplot_init.f90 \
	   math_routines.f90 \
	   geometryUVLM.f90 \
	   geometryModule.f90

OBJS_f901 = uvlm.o \
	   surfaceModule.o \
	   wakeModule.o \
	   solvers.o \
	   tecplot_init.o \
	   math_routines.o \
	   geometryUVLM.o \
	   geometryModule.o

SRCS_f902 = main.f90 \
	   ellipticModule.f90

OBJS_f902 = main.o \
	   ellipticModule.o
	   
SRC_DIR_f901 = ../source/
SRC_DIR_f902 =

OBJS_DIR = obj/
EXE_DIR = bin/

EXE = esfera
FC = ifort
IDIR = 
CFLAGS = -warn all -O2 -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include  -module $(OBJS_DIR) $(IDIR)

LFLAGS = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -lstdc++ 
LIBS = $(SRC_DIR_f901)tecio64.a 

VPATH = $(SRC_DIR_f901):$(SRC_DIR_f902):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f901) $(OBJS_f902))

all : $(EXE)
	@echo "All up to date!"

$(EXE) : $(OBJS_f901) $(OBJS_f902)
	@mkdir -p $(EXE_DIR)
	$(FC) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f901):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f901)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f902):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f902)$(@:.o=.f90) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)
	rm -f *.txt *.DAT *.dat *.plt

# Dependencies of files
uvlm.o: \
    uvlm.f90 \
    math_routines.o \
    surfaceModule.o
main.o: \
    main.f90 \
    math_routines.o \
    tecplot_init.o \
    solvers.o \
    geometryUVLM.o \
    ellipticModule.o \
    uvlm.o
tecplot_init.o: \
    tecplot_init.f90 \
    tecio.f90
math_routines.o: \
    math_routines.f90
surfaceModule.o: \
    surfaceModule.f90 \
    wakeModule.o \
    math_routines.o
wakeModule.o: \
    wakeModule.f90
solvers.o: \
    solvers.f90
geometryModule.o: \
	geometryModule.f90 \
	surfaceModule.o 
geometryUVLM.o: \
	geometryUVLM.f90 \
	uvlm.o \
	geometryModule.o \
	tecplot_init.o 
ellipticModule.o: \
	ellipticModule.f90 \
	surfaceModule.o
