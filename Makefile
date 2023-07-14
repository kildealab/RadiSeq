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
# -O Obtimize compilation
# -g Include debugging information
# -std=c++11 C++ language standard should be C++11
# -I$(INC_DIR) Include files from INC_DIR when compiling
CFLAGS = -c -Wall -O -g -std=c++11 -I$(INC_DIR)
LDFLAGS = 

SRC = $(wildcard $(SRC_DIR)/*.cpp)
# Take all source files from SRC_DIR
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC))
# For every source file, generate an object file in OBJ_DIR
EXECUTABLE = SDD_NGSSimulator
# Specifies the name of the executable generated


all: $(SRC) $(EXECUTABLE)
$(EXECUTABLE): $(OBJS)
	$(CPP) $(LDFLAGS) $(OBJS) -o $@
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) $(CFLAGS) $< -o $@
# $(SRC) $(EXECUTABLE) -- Dependencies for building the executable
# $(EXECUTABLE): $(OBJS) -- Rule for linking the object files
# rest are the commands: to link the object files and create the executable
# and to compile source files to object files respectively


# Remove all the object files and the executable
clean:
	rm -rf $(OBJ_DIR)/*.o $(EXECUTABLE)    