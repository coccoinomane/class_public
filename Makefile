#Some Makefile for CLASS.
#Julien Lesgourgues, 28.11.2011

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build

.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.c source:tools:main:test
vpath %.o build
vpath .base build

########################################################
###### LINES TO ADAPT TO YOUR PLATEFORM ################
########################################################

# your C compiler:
CC       = gcc
#CC       = icc
#CC       = pgcc

# your tool for creating static libraries:
AR        = ar rv

# (OPT) your python interpreter
PYTHON = python

# your optimization flag
OPTFLAG = -O4 -ffast-math #-march=native
#OPTFLAG = -Ofast -ffast-math #-march=native
#OPTFLAG = -fast

# your openmp flag (comment for compiling without openmp)
OMPFLAG   = -fopenmp
#OMPFLAG   = -mp -mp=nonuma -mp=allcores -g
#OMPFLAG   = -openmp

# all other compilation flags
CCFLAG = -g -fPIC -std=c99
#CFLAGS += -w
CCFLAG += -DDEBUG
LDFLAG = -g -fPIC

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. hyrec or ../hyrec)
HYREC = hyrec


########################################################
###### IN PRINCIPLE THE REST SHOULD BE LEFT UNCHANGED ##
########################################################

# pass current working directory to the code
CCFLAG += -D__CLASSDIR__='"$(MDIR)"'

# where to find include files *.h
INCLUDES = -I../include

# automatically add external programs if needed. First, initialize to blank.
EXTERNAL =

# directory where to put the documentation; a doc.html file will be created there
DOC_DIR = $(MDIR)/doc

# configuration file to generate the documentation (using doxygen)
DOXYGEN_FILE = $(DOC_DIR)/doxygen.txt

# Try to automatically avoid an error 'error: can't combine user with ...'
# which sometimes happens with brewed Python on OSX:
CFGFILE=$(shell $(PYTHON) -c "import sys; print sys.prefix+'/lib/'+'python'+'.'.join(['%i' % e for e in sys.version_info[0:2]])+'/distutils/distutils.cfg'")
PYTHONPREFIX=$(shell grep -s "prefix" $(CFGFILE))
ifeq ($(PYTHONPREFIX),)
PYTHONFLAGS=--user
else
PYTHONFLAGS=
endif

# eventually update flags for including HyRec
ifneq ($(HYREC),)
vpath %.c $(HYREC)
CCFLAG += -DHYREC
#LDFLAGS += -DHYREC
INCLUDES += -I../hyrec
EXTERNAL += hyrectools.o helium.o hydrogen.o history.o
endif

# If SONG is required, define the preprocessor macro WITH_SONG_SUPPORT
ifeq ($(WITH_SONG_SUPPORT),1)
	WITH_BISPECTRA = 1
	CCFLAG += -DWITH_SONG_SUPPORT
endif

# If bispectra or Fisher matrix are required, define the macro WITH_BISPECTRA
ifeq ($(WITH_BISPECTRA),1)
	CCFLAG += -DWITH_BISPECTRA
endif


%.o:  %.c .base
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o arrays.o parser.o quadrature.o hyperspherical.o common.o

SOURCE = input.o background.o thermodynamics.o perturbations.o primordial.o nonlinear.o transfer.o spectra.o lensing.o

INPUT = input.o

PRECISION = precision.o

BACKGROUND = background.o

THERMO = thermodynamics.o

PERTURBATIONS = perturbations.o

TRANSFER = transfer.o

PRIMORDIAL = primordial.o

SPECTRA = spectra.o

NONLINEAR = nonlinear.o

LENSING = lensing.o

OUTPUT = output.o

CLASS = class.o

TEST_LOOPS = test_loops.o

TEST_DEGENERACY = test_degeneracy.o

TEST_TRANSFER = test_transfer.o

TEST_NONLINEAR = test_nonlinear.o

TEST_PERTURBATIONS = test_perturbations.o

TEST_THERMODYNAMICS = test_thermodynamics.o

TEST_BACKGROUND = test_background.o

TEST_SIGMA = test_sigma.o

TEST_HYPERSPHERICAL = test_hyperspherical.o

TEST_STEPHANE = test_stephane.o

# Compile source files specific to SONG
ifeq ($(WITH_BISPECTRA),1)
	SOURCE += bessel.o bispectra.o fisher.o
	TOOLS += utility.o song_tools.o mesh_interpolation.o slatec_3j_C.o binary.o
endif

C_TOOLS =  $(addprefix tools/, $(addsuffix .c,$(basename $(TOOLS))))
C_SOURCE = $(addprefix source/, $(addsuffix .c,$(basename $(SOURCE) $(OUTPUT))))
C_TEST = $(addprefix test/, $(addsuffix .c,$(basename $(TEST_DEGENERACY) $(TEST_LOOPS) $(TEST_TRANSFER) $(TEST_NONLINEAR) $(TEST_PERTURBATIONS) $(TEST_THERMODYNAMICS))))
C_MAIN = $(addprefix main/, $(addsuffix .c,$(basename $(CLASS))))
C_ALL = $(C_MAIN) $(C_TOOLS) $(C_SOURCE)
H_ALL = $(addprefix include/, common.h svnversion.h $(addsuffix .h, $(basename $(notdir $(C_ALL)))))
PRE_ALL = cl_ref.pre clt_permille.pre
INI_ALL = explanatory.ini lcdm.ini
MISC_FILES = Makefile CPU psd_FD_single.dat myselection.dat myevolution.dat README bbn/sBBN.dat external_Pk/* cpp
PYTHON_FILES = python/classy.pyx python/setup.py python/cclassy.pxd python/test_class.py




all: class libclass.a classy

libclass.a: $(TOOLS) $(SOURCE) $(EXTERNAL)
	$(AR)  $@ $(addprefix build/, $(TOOLS) $(SOURCE) $(EXTERNAL))

class: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o class $(addprefix build/,$(notdir $^)) -lm

biclass: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o class $(addprefix build/,$(notdir $^)) -lm

test_sigma: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_SIGMA)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o test_sigma $(addprefix build/,$(notdir $^)) -lm

test_loops: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_stephane: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_STEPHANE)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_degeneracy: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_DEGENERACY)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_transfer: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_TRANSFER)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_nonlinear: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_NONLINEAR)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_perturbations: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_PERTURBATIONS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_thermodynamics: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_THERMODYNAMICS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_background: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_BACKGROUND)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_hyperspherical: $(TOOLS) $(TEST_HYPERSPHERICAL)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o test_hyperspherical $(addprefix build/,$(notdir $^)) -lm


tar: $(C_ALL) $(C_TEST) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)
	tar czvf class.tar.gz $(C_ALL) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)

classy: libclass.a python/classy.pyx python/cclassy.pxd
	cd python; export CC=$(CC); $(PYTHON) setup.py install $(PYTHONFLAGS)

clean: .base
	rm -rf $(WRKDIR);
	rm -f libclass.a
	rm -f $(MDIR)/python/classy.c
	rm -rf $(MDIR)/python/build

doxy:
	doxygen $(DOXYGEN_FILE)
	rm -rf $(DOC_DIR)/doc.html
	ln -s $(DOC_DIR)/doxygen/html/index.html $(DOC_DIR)/doc.html

cleandoxy:
	rm -rf $(DOXYGEN_DIR)
