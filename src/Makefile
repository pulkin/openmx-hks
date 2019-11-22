IDIR =../include
CC=gcc -g -Wall
CFLAGS+=-I$(IDIR)

ODIR=obj
LDIR=../lib
BDIR=../build

LIBS=-lm -lhdf5

_DEPS = openmx-hks-lib.h simplemat.h simplejson.h simpleh5.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = openmx-hks-lib.o simplemat.o openmx-hks.o simplejson.o simpleh5.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

openmx-hks: $(OBJ)
	gcc -g -o $(BDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ $(BDIR)/*
