C = gcc
OPTS = -c -fPIC
LOADER = gcc

SOURCES_C = ks.c 
DEPENDS = $(SOURCES_C:.c=.d)
OBJECTS = $(SOURCES_C:.c=.o) 

BoutrosLab.statistics.general: $(OBJECTS)
	$(R_HOME)/bin/R CMD SHLIB -o BoutrosLab.statistics.general.so $(OBJECTS)

clean:
	rm -f *.o *.so
