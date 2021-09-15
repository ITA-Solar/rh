# -----------------------------------------------------------------------
#
#   Makefile for RH main library
#
#   Users should not need to edit anything in this file.
#   Please look at Makefile.config to change compilers, flags, libraries
#
# -----------------------------------------------------------------------

include Makefile.config

.SUFFIXES: .f90

.f90.o:
	$(F90C) -c $(F90FLAGS) $<

DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE.c = $(CC) $(CFLAGS) $(DEPFLAGS) -D$(OS) -c
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

SRC = abundance.c \
	  accelerate.c \
	  backgropac_xdr.c \
	  background.c \
	  barklem.c \
	  broad.c \
	  brs_xdr.c \
	  chemequil.c \
	  cocollisions.c \
          collision.c \
	  complex.c \
	  cubeconvol.c \
	  duplicate.c \
	  error.c \
	  expint.c \
	  expspline.c \
	  fillgamma.c \
	  fixedrate.c \
	  fpehandler.c \
	  gammafunc.c \
	  gaussleg.c \
	  getcpu.c \
	  getlambda.c \
	  getline.c \
	  giigen.c \
	  h2collisions.c \
	  humlicek.c \
	  hunt.c \
	  hydrogen.c \
	  initial_xdr.c \
	  initscatter.c \
	  iterate.c \
	  kurucz.c \
	  linear.c \
	  ltepops.c \
	  ludcmp.c \
	  matrix.c \
	  maxchange.c \
	  metal.c \
	  molzeeman.c \
	  nemetals.c \
	  ohchbf.c \
	  opacity.c \
	  options.c \
	  order.c \
	  parse.c \
	  paschen.c \
	  planck.c \
	  pops_xdr.c \
	  profile.c \
	  radrate_xdr.c \
	  rayleigh.c \
	  readatom.c \
	  readb_xdr.c \
	  readinput.c \
	  readj.c \
	  readmolecule.c \
	  readvalue.c \
	  redistribute.c \
	  scatter.c \
	  solvene.c \
	  sortlambda.c \
	  spline.c \
	  statequil.c \
	  stokesopac.c \
	  stopreq.c \
	  thomson.c \
	  vacuumtoair.c \
	  voigt.c \
	  w3.c \
	  wigner.c \
	  writeatmos_xdr.c \
	  writeatom_xdr.c \
	  writecoll_xdr.c \
	  writedamp_xdr.c \
	  writeinput_xdr.c \
	  writemetal_xdr.c \
	  writemolec_xdr.c \
	  writeopac_xdr.c \
	  writespect_xdr.c \
	  zeeman.c \
          trilinear_interp.c \
          linspace.c

SRC_F = hui_.f90 humlicek_.f90

OBJS = $(SRC:.c=.o)
OBJS_F = $(SRC_F:.f90=.o)
DEPS = $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRC))))
LIBS = librh.a librh_f90.a

all: $(LIBS)

%.o : %.c
%.o : %.c $(DEPDIR)/%.d
	$(COMPILE.c) $<
	$(POSTCOMPILE)

librh.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@

librh_f90.a: $(OBJS_F)
	$(AR) $(ARFLAGS) $@ $(OBJS_F)
	$(RANLIB) $@

clean:
	rm -f $(OBJS) $(OBJS_F) $(LIBS)

-include $(DEPS)
