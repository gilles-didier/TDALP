# list of executable files to produce
TDALP = tdalp


# .o files necessary to build the executables
OBJ_TDALP = Utils.o StateTree.o Tree.o Fraction.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeGeneric.o DrawTreeResultGeneric.o ContinuousParsimony.o ContinuousReconstruction.o 
########### MODIFY ONLY TO CHANGE OPTIONS ############

# compiler and its options
CC = gcc
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -g -march=native -I/usr/include/cairo



# linker and its options
LD = $(CC)
LDFLAGS =  -lm  -lcairo -g


############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(DRAW) $(THEORETICAL_TEST) $(EMPIRICAL_TEST) $(EMPIRICAL_TEST_BIN) $(CONTROL_TEST)


# build the executable
$(TDALP): $(OBJ_TDALP)
	$(LD) $^ -o $@ $(LDFLAGS)

############ DO NOT MODIFY ANYTHING BELOW THIS LINE ##############

# create .o from .c
.c.o:
	$(CC) $(CFLAGS) -c $<

# remove non essential files
clean:
	$(RM) *.o *~ *.log Makefile.d

# clean everything but sources
distclean: clean
	$(RM) $(EXE)

# dependencies
Makefile.d:
	$(CC) -MM $(CFLAGS) *.c > Makefile.d

# only real files can be non phony targets
.PHONY: all clean distclean debug release

-include Makefile.d
