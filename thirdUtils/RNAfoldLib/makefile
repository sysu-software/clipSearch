CC=gcc
CXXC=g++
LIBS=-lm
INCLUDE=-I.
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE

O =  energy_par.o fold.o fold_vars.o utils.o params.o gquad.o aln_util.o

viennaRNA.a: $(O) 
	ar rcus libRNAfold.a $(O)
	ranlib libRNAfold.a
	
energy_par.o: energy_par.c
	$(CC) $(CFLAGS) $(HG_DEFS) -c energy_par.c

fold.o: fold.c fold.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c fold.c

fold_vars.o: fold_vars.c
	$(CC) $(CFLAGS) $(HG_DEFS) -c fold_vars.c

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c utils.c

gquad.o: gquad.c gquad.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c gquad.c

aln_util.o: aln_util.c aln_util.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c aln_util.c
	
params.o: params.c
	$(CC) $(CFLAGS) $(HG_DEFS) -c params.c
			
clean:
	rm -f *.o