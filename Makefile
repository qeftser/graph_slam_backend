
# compiler
CC := gcc

# compiler and linker flags
CC_FLAGS := -g
LINKER_FLAGS := -g -lm

# the source files for the program
SRC_FILES := src/graph_slam.c src/hash.c src/homogeneous.c src/solver.c src/set.c
HDR_FILES := src/graph_slam.h src/default.h src/homogeneous.h src/hash.h src/solver.h src/set.h
OBJ_FILES := $(SRC_FILES:.c=.o)

# prefix for install
PREFIX := /usr/local

# the library name 
ALIB := graph_slam_backend.a

all: $(OBJ_FILES)

test: $(OBJ_FILES) Makefile src/test.c
	mkdir -p bin && $(CC) $(OBJ_FILES) src/test.c $(CXX_FLAGS) $(LINKER_FLAGS) -o ./bin/$@.e && ./bin/test.e

scratch: $(OBJ_FILES) Makefile src/scratch.c
	mkdir -p bin && $(CC) $(OBJ_FILES) src/scratch.c $(CXX_FLAGS) $(LINKER_FLAGS) -o ./bin/$@.e

.c.o:
	$(CC) -c $(CXX_FLAGS) $< -o $@

clean:
	rm -f $(OBJ_FILES) ./bin/test.e ./bin/scratch.e ./bin/graph_slam_backend.a

# construct the archive file
archive: $(OBJ_FILES) 
	mkdir -p bin && ar rcs ./bin/$(ALIB) $(OBJ_FILES)

# assume a unix-like system
install: $(OBJ_FILES) archive
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	mkdir -p $(DESTDIR)$(PREFIX)/include/graph_slam_backend
	cp ./bin/$(ALIB) $(DESTDIR)$(PREFIX)/lib/$(ALIB)
	cp $(HDR_FILES) $(DESTDIR)$(PREFIX)/include/graph_slam_backend

# remove the files
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/lib/$(ALIB)
	rm -f -r $(DESTDIR)$(PREFIX)/include/graph_slam_backend
