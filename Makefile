# Makefile for Flight_Sim

# Deactivate implicit rules
.SUFFIXES:

# Directories
SRC_DIR = src
COM_DIR = common
BIN_DIR = bin

# List common files (ordered based on dependency)
COMMON_FILES = json.f90 jsonx.f90 linalg.f90 micro_time.f90 database_m.f90 udp_windows_m.f90 connection_m.f90
COMMON_PATHS = $(addprefix $(COM_DIR)/, $(COMMON_FILES))

# List source files (ordered based on dependency)
SRC_FILES = goates.f90 aircraft.f90 sim.f90
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
PROGRAM = sim.exe

# Default make
default:
	$(COMPILER) $(FLAGS) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH) -o $(PROGRAM)  -lws2_32

# Debug option
debug:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) $(DEBUG_FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH) 

# Debug with all warnings
wall:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) $(DEBUG_FLAGS) -Wall -fcheck=all -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH) 

# Cleanup
clean:
	rm -rf *.mod *.o *.dll *.exe $(SRC_DIR)/*.mod $(COM_DIR)/*.mod