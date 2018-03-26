# bacio Makefile
 
###############################################################
#
#   AUTHOR:    Gilbert - W/NP11
#
#   DATE:      01/11/1999
#
#   PURPOSE:   This script uses the make utility to update the bacio 
#              archive libraries.
#
###############################################################

include ../../../../macros.make

FCMP   = $(FCserial)
CCMP   = $(CC)
LIB    = ../../../../libbacio_4.a
FFLAGS = $(BACIO_FFLAGS)
AFLAGS = $(ARFLAGS)
CFLAGS = $(BACIO_CFLAGS)
OBJS   = bacio.o baciof.o bafrio.o byteswap.o chk_endianc.o

$(LIB): $(OBJS)
	ar -rv $(AFLAGS) $@ $^

bacio.o: bacio.c clib.h
	$(CCMP) -c $(CFLAGS) bacio.c

baciof.o: baciof.f
	$(FCMP) -c $(FFLAGS) baciof.f

bafrio.o: bafrio.f
	$(FCMP) -c $(FFLAGS) bafrio.f

byteswap.o: byteswap.c 
	$(CCMP) -c $(CFLAGS) byteswap.c

chk_endianc.o: chk_endianc.f 
	$(FCMP) -c $(FFLAGS) chk_endianc.f

clean:
	$(RM) *.o *.mod $(LIB)
