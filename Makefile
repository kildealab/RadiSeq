CPP = g++
# Specifies the C++ compiler

INC_DIR = include
SRC_DIR = source
OBJ_DIR = objects
# Directory where header files (.h) are stored
# Directory where source files (.cpp) are stored
# Destination directory to store object files (.o)


# CFLAGS are compiler flags.
# -c Compile one souce code at a time and link them only after all objext files are generated
# -Wall Enable "warning all" 
# -O3 Obtimize compilation to the 3rd degree (highest level)
# -g Include debugging information
# -std=c++1z or c++17 C++ language standard should be C++17 or above in clang or gcc respectively
# -fopenmp Use openmp library for parallelization
# -I$(INC_DIR) Include files from INC_DIR when compiling
# -MMD -MP: emit a .d file per .o listing its header dependencies (pulled in below via -include),
# 	so editing a header rebuilds every .o that includes it, not just the .cpp that changed.
# 	Without this, mismatched .o files can silently link together with different header layouts.
# -fno-strict-aliasing: needed for this codebase's raw pointer aliasing over mmap'd buffers;
# 	without it, -O3 has been observed to miscompile InduceSeq's genome-processing path.

# Determine the C++17 standard flag based on the compiler

COMPILER := $(shell $(CPP) -dM -E - < /dev/null | grep __clang__)

ifneq ($(COMPILER),)
    CXXFLAGS = -c -Wall -O3 -fno-strict-aliasing -g -std=c++1z -fopenmp -MMD -MP -I$(INC_DIR)
else
    CXXFLAGS = -c -Wall -O3 -fno-strict-aliasing -g -std=c++17 -fopenmp -MMD -MP -I$(INC_DIR)
endif

LDFLAGS = -fopenmp -lz
# LDFLAGS are linker flags
# -fopenmp link openmp library for parallelization
# -lrt link against the real-time extensions library (sometimes needed when working with POSIX functions)

SRC = $(wildcard $(SRC_DIR)/*.cpp)
# Take all source files from SRC_DIR
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC))
# For every source file, generate an object file in OBJ_DIR
EXECUTABLE = RadiSeq
# Specifies the name of the executable generated
EXECUTABLE2 = RadiSeqProfiler
# Seperate executable for the profiler program

# Source file for the new program
SRC2 = radiSeqProfiler/radiSeqProfiler.cpp
OBJ2 = $(OBJ_DIR)/radiSeqProfiler.o

# Create the 'objects' directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR))

all: $(EXECUTABLE) $(EXECUTABLE2)
$(EXECUTABLE): $(OBJS)
	$(CPP) $(OBJS) -o $@ $(LDFLAGS)
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) $(CXXFLAGS) $< -o $@
# $(SRC) $(EXECUTABLE) -- Dependencies for building the executable
# $(EXECUTABLE): $(OBJS) -- Rule for linking the object files
# rest are the commands: to link the object files and create the executable
# and to compile source files to object files respectively

# Rule for building the new executable
$(EXECUTABLE2): $(OBJ2)
	$(CPP) $(OBJ2) -o $@ $(LDFLAGS)
# Rule for compiling the new radiSeqProfiler.cpp file
$(OBJ_DIR)/radiSeqProfiler.o: $(SRC2)
	$(CPP) $(CXXFLAGS) $< -o $@

# Pull in the auto-generated header dependency rules (see the -MMD -MP note above), so editing
# a header correctly triggers a rebuild of every .o that (transitively) includes it. Silently
# does nothing on a clean checkout / after `make clean`, before any .d files exist yet.
-include $(OBJS:.o=.d) $(OBJ2:.o=.d)

# Remove all the object files and the executable
clean:
	rm -rf $(OBJ_DIR)/*.o $(OBJ_DIR)/*.d $(EXECUTABLE) $(EXECUTABLE2)