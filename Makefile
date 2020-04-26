#Set to your system's CPLEX path, architecture, etc. (or make sure the appropriate paths are in CPATH and LIBRARY_PATH)
ILOGBASE = /opt/ibm/ILOG/CPLEX_Studio128
ILOGSYSTEM = x86-64_linux
ILOGLIBFORMAT = static_pic
ILOGLIBSUFFIX = $(ILOGSYSTEM)/$(ILOGLIBFORMAT)
#BLAS library name as installed
BLASLIB = openblaso

ILOGFLAGS = -DILOGLUE -DIL_STD -I$(ILOGBASE)/cpoptimizer/include/ -I$(ILOGBASE)/concert/include/ -I$(ILOGBASE)/cplex/include/ 
ILOGLIBS  = -L$(ILOGBASE)/cpoptimizer/lib/$(ILOGLIBSUFFIX) -L$(ILOGBASE)/concert/lib/$(ILOGLIBSUFFIX) -L$(ILOGBASE)/cplex/lib/$(ILOGLIBSUFFIX) -lcp -lilocplex -lcplex -lconcert 
ARFLAGS = -fopenmp
ARLIBS = -lm -larmadillo -l$(BLASLIB) -ldl

IDIR=.

CC=gcc
CXX=g++
CFLAGS=-g -Wall -w -O3 -I$(IDIR) $(ILOGFLAGS) $(ARFLAGS) -I $(IDIR)/include 

ODIR=.

LIBS=-L. $(ILOGLIBS) $(ARLIBS)

_DEPS = AgileCombiFD.h PostProcessor.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = AgileCombiFD.o PostProcessor.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all:     main

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

main: $(OBJ) $(ODIR)/main.o
	$(CXX) -o iafd $(OBJ) $(ODIR)/main.o $(CFLAGS) $(LIBS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ iafd
