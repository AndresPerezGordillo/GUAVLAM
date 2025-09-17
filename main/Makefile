#
# Makefile for code 'GUAVLAM' 
# Author: Juan D. Colmenares
#

SRCS_f901 = uvlm.f90 \
	   surfaceModule.f90 \
	   wakeModule.f90 \
	   tecplot_init.f90 \
	   math_routines.f90 \
	   uvlm_wrappers.f90 \
	   caseModule.f90

OBJS_f901 = uvlm.o \
	   surfaceModule.o \
	   wakeModule.o \
	   tecplot_init.o \
	   math_routines.o \
	   uvlm_wrappers.o \
	   caseModule.o

SRCS_f902 = main.f90 

OBJS_f902 = main.o 
	   
SRCS_f903 = rectangularBladeModule.f90 \
	   ellipticModule.f90 \
	   crearRotor.f90 

OBJS_f903 = rectangularBladeModule.o \
	   ellipticModule.o \
	   crearRotor.o 
	   
SRC_DIR_f901 = ../source/
SRC_DIR_f902 =
SRC_DIR_f903 =

OBJS_DIR = obj/
EXE_DIR = bin/
ACML = /home/andres/amdsdk/v1.0/acml/open64_64

EXE = guavlam
FC = openf95
IDIR = 
OMPFLAG = -mp
CFLAGS = -Wall $(OMPFLAG) -march=auto -module $(OBJS_DIR) $(IDIR)
LFLAGS = -lstdc++
LIBS = $(SRC_DIR_f901)tecio64.a $(ACML)/lib/libacml.a

VPATH = $(SRC_DIR_f901):$(SRC_DIR_f902):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f901) $(OBJS_f902))

all : $(EXE)
	@echo "All up to date!"

$(EXE) : $(OBJS_f901) $(OBJS_f902)
	@mkdir -p $(EXE_DIR)
	$(FC) -Wall $(OMPFLAG) -o $(EXE_DIR)$(EXE) $(OBJS) $(LIBS) $(LFLAGS)

$(OBJS_f901):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f901)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f902):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f902)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f903):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f903)$(@:.o=.f90) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)
	rm -f *.txt *.DAT *.plt tp1* *.geom *.wake CaseControl*

remake :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)
	@make all

clear :
	rm -f *.txt *.DAT *.plt tp1* *.geom *.wake CaseControl*

geom : PARAMETERS.PAR geom.exe
	./geom.exe

geom.exe : $(OBJS_f903)
	$(FC) $(OMPFLAG) -o geom.exe $(OBJS_DIR)math_routines.o $(OBJS_DIR)wakeModule.o \
		$(OBJS_DIR)surfaceModule.o $(OBJS_DIR)rectangularBladeModule.o \
		$(OBJS_DIR)ellipticModule.o $(OBJS_DIR)crearRotor.o $(LIBS) -ipa

# Dependencies of files
crearRotor.o : \
	surfaceModule.o \
	rectangularBladeModule.o \
	ellipticModule.o \
	crearRotor.f90

uvlm.o: \
    uvlm.f90 \
    math_routines.o \
    surfaceModule.o
main.o: \
    main.f90 \
    uvlm_wrappers.o 
tecplot_init.o: \
    tecplot_init.f90 \
    tecio.f90
math_routines.o: \
    math_routines.f90
rectangularBladeModule.o: \
    rectangularBladeModule.f90 \
    math_routines.o \
    surfaceModule.o
surfaceModule.o: \
    surfaceModule.f90 \
    wakeModule.o \
    math_routines.o
wakeModule.o: \
    wakeModule.f90
uvlm_wrappers.o: \
	uvlm_wrappers.f90 \
	uvlm.o \
	tecplot_init.o 
ellipticModule.o: \
	ellipticModule.f90 \
	surfaceModule.o
caseModule.o: \
	caseModule.f90 \
	surfaceModule.o 
