# -----------------------------------------------------------------------
#
#   Makefile for RH 1.5D
#
#   Users should not need to edit anything in this file.
#   Please look at Makefile.config to change compilers, flags, libraries
#
# -----------------------------------------------------------------------

include ../Makefile.config

# Try to get git revision, set preprocessor macros if successful
REV=$(shell git describe --always)
ifeq ($(REV), )
    CFLAGS += -I../ -I$(HDF5_DIR)/include/
else
    REV = $(shell git log -1 --pretty=format:"%h  %an  %ai")
    CFLAGS += -I../ -I$(HDF5_DIR)/include/ -D'REV_ID="$(REV)"'
endif

LDFLAGS += -L../ -L$(HDF5_DIR)/lib/ -lrh -lrh_f90 -lhdf5 -lhdf5_hl -lm -lpthread

SRC = accelerate_p.c  \
  	  anglequad.c  \
	  background_p.c  \
      bezier.c  \
	  distribute_jobs.c  \
	  escape.c  \
	  expint_p.c \
	  feautrier.c  \
	  formal.c  \
	  hdf5atmos.c  \
	  hydrostat.c  \
	  initial_p.c  \
	  iterate_p.c  \
	  ludcmp_p.c  \
	  multiatmos.c  \
	  parallel.c  \
	  piecestokes.c  \
	  piecewise.c  \
	  project.c  \
	  readRayInput.c  \
	  readatmos.c  \
	  redistribute_p.c  \
	  riiplane.c  \
	  scatter_p.c  \
	  statequil_p.c  \
	  writeAux_p.c  \
	  writeRay.c  \
	  writeindata_p.c \
          reademistab.c \
          trapezoidal.c \
          readpopsin.c \
          hdf5popsin.c

BIN_SRC = rh15d_ray.c rh15d_ray_pool.c rh15d_lteray.c
OBJS = $(SRC:.c=.o)
BIN_OBJ = $(BIN_SRC:.c=.o)

# Dependencies and build rule
DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
COMPILE.c = $(CC) $(DEPFLAGS) $(CFLAGS) -c
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

%.o : %.c
%.o : %.c $(DEPDIR)/%.d
	$(COMPILE.c) $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d
DEPS = $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(OBJS) $(BIN_OBJ))))

# Rules for the executables
BIN = $(basename $(BIN_OBJ))

all: $(BIN)

$(BIN): % : $(OBJS) %.o ../librh.a ../librh_f90.a
	$(LD) -o $@ $(OBJS) $@.o $(LDFLAGS)

../librh.a ../librh_f90.a:
	cd ..; $(MAKE) librh.a  librh_f90.a

clean:
	rm -f $(OBJS) $(BIN_OBJ) $(BIN)

-include $(DEPS)
