# =============================================================
# Makefile for the Unsteady DVM project
# Compiler: gfortran (change FC for ifort, etc.)
#
# Usage:
#   make           – compile everything
#   make run       – compile and run with config.nml
#   make clean     – remove build artifacts and executable
#   make veryclean – also remove output files
# =============================================================

FC      = gfortran
# -O2: optimise | -Wall -Wextra: warnings | -fcheck=all: runtime checks
# Remove -fcheck=all for production runs (slower but catches bugs during dev)
FFLAGS  = -O2 -Wall -Wextra -fcheck=all -g

SRCDIR  = src
OBJDIR  = obj
TARGET  = udvm

# Source files MUST be listed in dependency order (used modules before using modules)
SRCS =  parameters.f90  \
        kinematics.f90   \
        geometry.f90     \
        influence.f90    \
        wake.f90         \
        forces.f90       \
        io.f90           \
        solver.f90       \
        main.f90

OBJS = $(patsubst %.f90,$(OBJDIR)/%.o,$(SRCS))

# ---- Targets ------------------------------------------------

.PHONY: all clean veryclean run

all: $(OBJDIR) output $(TARGET)

# Create object and output directories if they don't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

output:
	mkdir -p output

# Link
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo ""
	@echo "  Build successful: ./$(TARGET)"
	@echo ""

# Compile each source file
# -J$(OBJDIR)  : put generated .mod files in obj/
# -I$(OBJDIR)  : look for .mod files in obj/
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

# ---- Explicit dependencies (for make -j safety) -------------
$(OBJDIR)/kinematics.o: $(OBJDIR)/parameters.o
$(OBJDIR)/geometry.o:   $(OBJDIR)/parameters.o
$(OBJDIR)/influence.o:  $(OBJDIR)/parameters.o
$(OBJDIR)/wake.o:       $(OBJDIR)/parameters.o
$(OBJDIR)/forces.o:     $(OBJDIR)/parameters.o
$(OBJDIR)/io.o:         $(OBJDIR)/parameters.o
$(OBJDIR)/solver.o:     $(OBJDIR)/parameters.o $(OBJDIR)/influence.o $(OBJDIR)/wake.o
$(OBJDIR)/main.o:       $(OBJDIR)/parameters.o $(OBJDIR)/geometry.o  \
                        $(OBJDIR)/kinematics.o  $(OBJDIR)/influence.o \
                        $(OBJDIR)/wake.o        $(OBJDIR)/solver.o    \
                        $(OBJDIR)/forces.o      $(OBJDIR)/io.o

# ---- Convenience targets ------------------------------------

run: all
	./$(TARGET) config.nml

clean:
	rm -rf $(OBJDIR) $(TARGET)

veryclean: clean
	rm -rf output