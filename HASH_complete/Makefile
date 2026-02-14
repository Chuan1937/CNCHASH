# Makefile for HASH v1.2
# Updated for new directory structure

# Fortran compiler (change to your preference)
FC = gfortran
FFLAGS = -O

# Directories
SRCDIR = src
SUBDIR = $(SRCDIR)/subs
DRVDIR = $(SRCDIR)/drivers
INCDIR = $(SRCDIR)/include
OBJDIR = obj
BINDIR = bin

# Include path
INCS = -I$(INCDIR)

# Object files for each driver
OBJS1 = $(OBJDIR)/fmech_subs.o \
        $(OBJDIR)/pol_subs.o \
        $(OBJDIR)/station_subs.o \
        $(OBJDIR)/uncert_subs.o \
        $(OBJDIR)/util_subs.o

OBJS2 = $(OBJDIR)/fmech_subs.o \
        $(OBJDIR)/pol_subs.o \
        $(OBJDIR)/station_subs.o \
        $(OBJDIR)/uncert_subs.o \
        $(OBJDIR)/util_subs.o \
        $(OBJDIR)/vel_subs.o

OBJS3 = $(OBJDIR)/fmamp_subs.o \
        $(OBJDIR)/pol_subs.o \
        $(OBJDIR)/station_subs.o \
        $(OBJDIR)/uncert_subs.o \
        $(OBJDIR)/util_subs.o \
        $(OBJDIR)/vel_subs.o

OBJS4 = $(OBJDIR)/fmech_subs.o \
        $(OBJDIR)/pol_subs.o \
        $(OBJDIR)/station_subs_5char.o \
        $(OBJDIR)/uncert_subs.o \
        $(OBJDIR)/util_subs.o \
        $(OBJDIR)/vel_subs.o

# Default target
all: dirs $(BINDIR)/hash_driver1 $(BINDIR)/hash_driver2 $(BINDIR)/hash_driver3 $(BINDIR)/hash_driver4 $(BINDIR)/hash_driver5

dirs:
	mkdir -p $(OBJDIR) $(BINDIR)

# Drivers
$(BINDIR)/hash_driver1: $(DRVDIR)/hash_driver1.f $(OBJS1)
	$(FC) $(FFLAGS) $(INCS) $< $(OBJS1) -o $@

$(BINDIR)/hash_driver2: $(DRVDIR)/hash_driver2.f $(OBJS2)
	$(FC) $(FFLAGS) $(INCS) $< $(OBJS2) -o $@

$(BINDIR)/hash_driver3: $(DRVDIR)/hash_driver3.f $(OBJS3)
	$(FC) $(FFLAGS) $(INCS) $< $(OBJS3) -o $@

$(BINDIR)/hash_driver4: $(DRVDIR)/hash_driver4.f $(OBJS4)
	$(FC) $(FFLAGS) $(INCS) $< $(OBJS4) -o $@

$(BINDIR)/hash_driver5: $(DRVDIR)/hash_driver5.f $(OBJS2)
	$(FC) $(FFLAGS) $(INCS) $< $(OBJS2) -o $@

# Subroutines
$(OBJDIR)/pol_subs.o: $(SUBDIR)/pol_subs.f
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(OBJDIR)/fmech_subs.o: $(SUBDIR)/fmech_subs.f
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(OBJDIR)/fmamp_subs.o: $(SUBDIR)/fmamp_subs.f
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(OBJDIR)/station_subs.o: $(SUBDIR)/station_subs.f
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(OBJDIR)/station_subs_5char.o: $(SUBDIR)/station_subs_5char.f
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(OBJDIR)/uncert_subs.o: $(SUBDIR)/uncert_subs.f
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(OBJDIR)/util_subs.o: $(SUBDIR)/util_subs.f
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(OBJDIR)/vel_subs.o: $(SUBDIR)/vel_subs.f
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

# Clean
clean:
	rm -rf $(OBJDIR) $(BINDIR)

# Install (optional)
install: all
	cp -r data examples $(BINDIR)/

.PHONY: all dirs clean install
