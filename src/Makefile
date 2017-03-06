IDIR =../include
CC=gcc -g
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = openmx-hks-lib.h simplemat.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = openmx-hks-lib.o simplemat.o openmx-hks.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

openmx-hks: $(OBJ)
	gcc -g -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
