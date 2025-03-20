########################################################################
####################### Makefile Template ##############################
########################################################################
# Compiler settings - Can be customized.
CC = mpif90
CPP = mpif90 -cpp
# -I../MODULES/ tells the compiler where to look for .mod files,
# -J../MODULES tells it where to place generated .mod files.
CXXFLAGS = -g -O0 -ffree-line-length-0 -march=native -funroll-loops -fopenmp -I../MODULES/ -J../MODULES
LDFLAGS = -lblas -llapack

# Makefile settings - Can be customized.
APPNAME = myapp
EXT = .f90
SRCDIR = src
OBJDIR = obj
MODULESDIR = ../MODULES

############## Do not change anything from here downwards! #############
SRC = $(wildcard $(SRCDIR)/*$(EXT))
OBJ = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)/%.o)
DEP = $(OBJ:$(OBJDIR)/%.o=%.d)

# Module source files and their object files.
MODULE_SRCS = $(wildcard $(MODULESDIR)/*$(EXT))
MODULE_OBJS = $(MODULE_SRCS:$(MODULESDIR)/%$(EXT)=$(OBJDIR)/%.o)

# UNIX-based OS variables & settings
RM = rm
DELOBJ = $(OBJ) $(MODULE_OBJS)
# Windows OS variables & settings
DEL = del
EXE = .exe
WDELOBJ = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)\\%.o)

########################################################################
####################### Targets beginning here #########################
########################################################################

all: $(APPNAME)

# Build the final executable by linking both main and module objects.
$(APPNAME): $(OBJ) $(MODULE_OBJS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Ensure that the main object depends on the module objects.
# This forces the modules to be compiled (and their .mod files generated)
# before the main source file is compiled.
$(OBJDIR)/b-field-weyl_points.o: $(MODULE_OBJS)

# Create dependency rules for src files.
# %.d: $(SRCDIR)/%$(EXT)
# 	@$(CPP) $(CXXLAGS) $< -MM -MT $(@:%.d=$(OBJDIR)/%.o) > $@

# Include dependency files.
# -include $(DEP)

# Rule to build .o files from source files in src.
$(OBJDIR)/%.o: $(SRCDIR)/%$(EXT)
	@mkdir -p $(OBJDIR)
	$(CC) $(CXXFLAGS) -c $< -o $@

# Rule to build .o files from module source files in MODULESDIR.
$(OBJDIR)/%.o: $(MODULESDIR)/%$(EXT)
	@mkdir -p $(OBJDIR)
	$(CC) $(CXXFLAGS) -c $< -o $@

################### Additional target for running MPI program ##############
.PHONY: run
run: $(APPNAME)
	mpirun -np 2 ./$(APPNAME)

################### Cleaning rules for Unix-based OS ###################
.PHONY: clean
clean:
	$(RM) $(DELOBJ) $(DEP) $(APPNAME)

.PHONY: cleandep
cleandep:
	$(RM) $(DEP)

#################### Cleaning rules for Windows OS #####################
.PHONY: cleanw
cleanw:
	$(DEL) $(WDELOBJ) $(DEP) $(APPNAME)$(EXE)

.PHONY: cleandepw
cleandepw:
	$(DEL) $(DEP)
