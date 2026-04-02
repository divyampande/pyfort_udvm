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

# --- OS Detection ---
ifeq ($(OS),Windows_NT)
    MKDIR_OBJ = if not exist $(OBJDIR) mkdir $(OBJDIR)
    MKDIR_OUT = if not exist output mkdir output
    RM_DIR    = if exist $(OBJDIR) rmdir /S /Q $(OBJDIR)
    RM_OUT    = if exist output rmdir /S /Q output
    RM_OUTSP  = if exist output_sweep rmdir /S /Q output_sweep
    RM_EXE    = if exist $(TARGET).exe del /Q /F $(TARGET).exe
    EXEC      = $(TARGET).exe
    RUN_CMD   = $(EXEC) config.nml
else
    MKDIR_OBJ = mkdir -p $(OBJDIR)
    MKDIR_OUT = mkdir -p output
    RM_DIR    = rm -rf $(OBJDIR)
    RM_OUT    = rm -rf output
    RM_EXE    = rm -f $(TARGET)
    EXEC      = $(TARGET)
    RUN_CMD   = ./$(EXEC) config.nml
endif

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
	$(MKDIR_OBJ)

output:
	$(MKDIR_OUT)

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
	$(RUN_CMD)

clean:
	$(RM_DIR)
	$(RM_EXE)

veryclean: clean
	$(RM_OUT)
	$(RM_OUTSP)