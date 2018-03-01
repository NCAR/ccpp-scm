# libsrc/w3nco makefile template
#
# Build Double Precision (Size of Real 8-byte and default Integer) version 
# of libw3nco_d.a
#
#
include ../../../../macros.make
 
FC     = $(FCserial)
LIB    = ../../../../libw3nco_d.a
FFLAGS = $(W3NCO_FFLAGS)
AFLAGS = $(W3NCO_ARFLAGS)
CFLAGS = $(W3NCO_CFLAGS)
SRCS   = $(wildcard *.f *.c)
OBJS   = $(addsuffix .o, $(basename $(SRCS)))

%.o: %.f
	$(FC) -c $(FFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) $<

$(LIB): $(OBJS)
	ar $(AFLAGS) $@ $^

clean:
	$(RM) *.o *.mod $(LIB)
