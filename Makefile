# Copyright 2009, 2010, 2014 IPOL Image Processing On Line http://www.ipol.im/
# Author: Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
# Author: Bruno Galerne <bruno.galerne@parisdescartes.fr>
# Author: Lara Raad <lara.raad@cmla.ens-cachan.fr>

# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the copyright notice and this notice are preserved.  This file is
# offered as-is, without any warranty.

SRC	= quilting_main.c io_png.c quilting.c common.c patch_search.c \
          blending.c boundary_cut.c random_number.c
BIN	= quilting
COPT	= -std=c99 -O3
OPENMP  =

default: $(BIN)

openmp: OPENMP=-fopenmp
openmp: default

$(BIN): $(SRC)
	$(CC) $(COPT) $(OPENMP) -o $@ $^ -lpng -lfftw3l -lm

clean:
	$(RM) $(BIN)
