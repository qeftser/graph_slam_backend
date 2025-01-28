
# compiler
CC := gcc

# compiler and linker flags
CC_FLAGS := -g
LINKER_FLAGS := -g -lm

# the source files for the program
SRC_FILES := src/graph_slam.c src/hash.c src/homogeneous.c src/solver.c
OBJ_FILES := $(SRC_FILES:.c=.o)

all: test scratch

test: $(OBJ_FILES) Makefile src/test.c
	$(CC) $(OBJ_FILES) src/test.c $(CXX_FLAGS) $(LINKER_FLAGS) -o ./bin/$@.e

scratch: $(OBJ_FILES) Makefile src/scratch.c
	$(CC) $(OBJ_FILES) src/scratch.c $(CXX_FLAGS) $(LINKER_FLAGS) -o ./bin/$@.e

.c.o:
	$(CC) -c $(CXX_FLAGS) $< -o $@

clean:
	rm -f $(OBJ_FILES) ./bin/test.e ./bin/scratch.e
