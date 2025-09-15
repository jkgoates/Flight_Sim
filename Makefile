
# make for HI-Mach. A Modified-Newtonian Panel Method 

# Deactivate implicit rules
.SUFFIXES:

# Directories
SRC_DIR = src
COM_DIR = common
BIN_DIR = bin

# List common files (ordered based on dependency)
COMMON_FILES = helpers.f95 linked_list.f90 json.f95 json_xtnsn.f95 linalg.f90 math.f90 std_atmosphere.f90 timer.f95 quartic.f90 props_reader.f90 sort.f90
COMMON_OUT = helpers.o linked_list.o json.o json_xtnsn.o linalg.o math.o std_atmosphere.o timer.o quartic.o props_reader.o sort.o
COMMON_PATHS = $(addprefix $(COM_DIR)/, $(COMMON_FILES))

# List source files (ordered based on dependency)
SRC_FILES = flow.f90 base_geom.f90 panel.f90 stl.f90 vtk.f90 tri.f90 surface_mesh.f90 aero_solver.f90 heat_solver.f90 trajectory_handler.f90 himach.f90
SRC_OUT = flow.o base_geom.o panel.o stl.o vtk.o tri.o surface_mesh.o aero_solver.o heat_solver.o trajectory_handler.o himach.o
SRC_PATHS = $(addprefix $(SRC_DIR)/, $(SRC_FILES))

# Main
MAIN_PATH = src/main.f90

# Compiler
COMPILER = gfortran

# Flags
DEBUG_FLAGS = -fbounds-check -fbacktrace -g
OMP_FLAG = -fopenmp
FLAGS = -O2 -fdefault-real-8

# Program name
PROGRAM = himach.exe

# Default make
default:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH)

# Debug option
debug:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) $(DEBUG_FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH) 

# Debug with all warnings
wall:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) $(DEBUG_FLAGS) -Wall -fcheck=all -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH) 

# Cleanup
clean:
	rm -rf *.mod *.o *.dll *.exe $(SRC_DIR)/*.mod $(COM_DIR)/*.mod