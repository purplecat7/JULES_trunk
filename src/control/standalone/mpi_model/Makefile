#this is not automatically generated!
      
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(mpi_comms_mod.o)

all : $(OBJS) 

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(mpi_comms_mod.o): $(ARC)(logging_mod.o) $(ARC)(prognostics.o) $(ARC)(string_utils_mod.o)

# End of file.
