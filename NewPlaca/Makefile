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
	   uvlm_wrappers.f90 \
	   caseModule.f90

OBJS_f901 = uvlm.o \
	   surfaceModule.o \
	   wakeModule.o \
	   solvers.o \
	   tecplot_init.o \
	   math_routines.o \
	   uvlm_wrappers.o \
	   caseModule.o

SRCS_f902 = main.f90 \
	   rectangularWingModule.f90 

OBJS_f902 = main.o \
	   rectangularWingModule.o 
	   
SRC_DIR_f901 = ../source/
SRC_DIR_f902 =

OBJS_DIR = obj/
EXE_DIR = bin/

EXE = placa_plana
FC = openf95
IDIR = 
CFLAGS = -mp -O2 -module $(OBJS_DIR) $(IDIR)
LFLAGS = -lstdc++ 
LIBS = $(SRC_DIR_f901)tecio64.a $(ACML)/lib/libacml.a

VPATH = $(SRC_DIR_f901):$(SRC_DIR_f902):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f901) $(OBJS_f902))

all : $(EXE)
	@echo "All up to date!"

$(EXE) : geom $(OBJS_f901) $(OBJS_f902)
	@mkdir -p $(EXE_DIR)
	$(FC) $(CFLAGS) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f901):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f901)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f902):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f902)$(@:.o=.f90) -o $(OBJS_DIR)$@

geom : geom.exe CaseControl
	./geom.exe
	
geom.exe : math_routines.o wakeModule.o \
		surfaceModule.o rectangularWingModule.o crearPlaca.f90
	$(FC) $(CFLAGS) -c crearPlaca.f90 -o $(OBJS_DIR)crearPlaca.o
	$(FC) -mp -o geom.exe $(OBJS_DIR)math_routines.o $(OBJS_DIR)wakeModule.o \
		$(OBJS_DIR)surfaceModule.o $(OBJS_DIR)rectangularWingModule.o $(OBJS_DIR)crearPlaca.o \
		$(LIBS)

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)
	rm -f *.txt *.DAT *.plt *.geom *.wake *.exe CaseControl_*

clear :
	rm -f *.txt *.DAT *.plt *.geom *.wake tp1* CaseControl_*

# Dependencies of files

uvlm.o: \
    uvlm.f90 \
    math_routines.o \
    surfaceModule.o
main.o: \
    main.f90 \
    solvers.o \
    uvlm_wrappers.o \
    caseModule.o 
tecplot_init.o: \
    tecplot_init.f90 \
    tecio.f90
math_routines.o: \
    math_routines.f90
rectangularWingModule.o: \
    rectangularWingModule.f90 \
    math_routines.o \
    surfaceModule.o
surfaceModule.o: \
    surfaceModule.f90 \
    wakeModule.o \
    math_routines.o
wakeModule.o: \
    wakeModule.f90
solvers.o: \
    solvers.f90 \
    math_routines.o
uvlm_wrappers.o: \
	uvlm_wrappers.f90 \
	tecplot_init.o \
	uvlm.o
caseModule.o: \
	caseModule.f90 \
	surfaceModule.o 
