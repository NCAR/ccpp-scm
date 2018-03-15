#sp Makefile

###############################################################
#
#   AUTHOR:    Vuong - W/NP11
#
#   DATE:      12/04/2000
#
#   PURPOSE:   This script uses the make utility to update the libsp 
#              archive libraries.
#              It first reads a list of source files in the library and
#              then generates a makefile used to update the archive
#              libraries.  The make command is then executed for each
#              archive library, where the archive library name and 
#              compilation flags are passed to the makefile through 
#              environment variables.
#
#   REMARKS:   Only source files that have been modified since the last
#              library update are recompiled and replaced in the object
#              archive libraries.  The make utility determines this
#              from the file modification times.
#
#              New source files are also compiled and added to the object 
#              archive libraries
#
# 11-29-12     Mirvis: Added  Intel comp, optimization and OMP flags 
#		and libsp_*g.a builds
#
# 05-29-14     Mirvis: Updated to ver 2.0.2
###############################################################

#
#     Generate a list of object files that corresponds to the
#     list of Fortran ( .f ) files in the current directory
 
#     Set flags for Double Precision (Size of Real 8-byte and default Integer) version
#     of libsp_d.a
 
include ../../../../macros.make

FC     = $(FCserial)
SRCS   = $(wildcard *.f *.F)
OBJS   = $(addsuffix .o, $(basename $(SRCS)))
LIB    = ../../../../libsp_v2.0.2_d.a
# libsp_v2.0.2_d.a
FFLAGS = $(SP_FFLAGS)
AFLAGS = $(SP_ARFLAGS)

$(LIB): $(OBJS)
	$(AR) $(AFLAGS) $@ $^

.f.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(AFLAGS)  $@ $*.o
	$(RM) $*.o

.F.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(AFLAGS)  $@ $*.o
	$(RM) $*.o

clean:
	$(RM) *.o *.mod $(LIB)
