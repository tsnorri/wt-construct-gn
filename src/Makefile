include ../local.mk
include ../common.mk

OBJECTS = cmdline.o file_handling.o main.o 

all: construct-wt
	
clean:
	$(RM) $(OBJECTS) cmdline.c cmdline.h construct-wt
	
main.c: cmdline.c

construct-wt: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../lib/sdsl/build/lib/libsdsl.a $(LIBDISPATCH_LIBS) $(BOOST_LIBS)
