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

# Determine the C++17 standard flag based on the compiler
COMPILER := $(shell $(CPP) -dM -E - < /dev/null | grep __clang__)
ifneq ($(COMPILER),)
    CXXFLAGS = -c -Wall -O3 -g -std=c++1z -fopenmp -I$(INC_DIR)
else
    CXXFLAGS = -c -Wall -O3 -g -std=c++17 -fopenmp -I$(INC_DIR)
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

# Create the 'objects' directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR))

all: $(SRC) $(EXECUTABLE)
$(EXECUTABLE): $(OBJS)
	$(CPP) $(LDFLAGS) $(OBJS) -o $@
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) $(CXXFLAGS) $< -o $@
# $(SRC) $(EXECUTABLE) -- Dependencies for building the executable
# $(EXECUTABLE): $(OBJS) -- Rule for linking the object files
# rest are the commands: to link the object files and create the executable
# and to compile source files to object files respectively


# Remove all the object files and the executable
clean:
	rm -rf $(OBJ_DIR)/*.o $(EXECUTABLE)    