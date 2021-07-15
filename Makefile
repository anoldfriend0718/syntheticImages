SRC_DIR = .
VPATH = .
BUILD_DIR = .
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE


all:  library progs 
qxu:  synthesizeCokedStructure2D generateCases3D createMask2D createMask3D

tmp: 
	g++ -g -std=c++11 tmp.c -o $@ || exit 1

directories:
	if [ ! -d "include" ]; then mkdir include; fi
	if [ ! -d "lib" ]; then mkdir lib; fi

includes:  directories
	cp -f config/QSSLIB_config.h include
	cd src; make $@ || exit 1


library:  directories includes 
	cd src; make $@ || exit 1
	make lsmpqs_lib || exit 1

lsmpqs_lib:
	find src -name "*.o" > objs_file.tmp
	cat objs_file.tmp | xargs ar -ru $(BUILD_DIR)/lib/liblsmpqs.a 
	ranlib lib/liblsmpqs.a
	rm -f objs_file.tmp 

progs:
	cd src; make $@ || exit 1

synthesizeCokedStructure2D:
	gcc $(DEBUG_FLAGS) $(QSSLIB_DIRS) \
	-I$(QSSLIB_INCLUDE) ./src/applications/synthesizeCokedStructure2D.c \
	-o $@ $(QSSLIB_LIBS) $(OMP_FLAGS)
	mv -f $@ $(QSSLIB_BIN_DIR) || exit 1

generateCases3D:
	gcc $(DEBUG_FLAGS) $(QSSLIB_DIRS) \
	-I$(QSSLIB_INCLUDE) ./src/applications/generateCases3D.c \
	-o $@ $(QSSLIB_LIBS) $(OMP_FLAGS)
	mv -f $@ $(QSSLIB_BIN_DIR) || exit 1

createMask2D:
	gcc $(DEBUG_FLAGS) $(QSSLIB_DIRS) \
	-I$(QSSLIB_INCLUDE) ./src/applications/createMask2D.c \
	-o $@ $(QSSLIB_LIBS) $(OMP_FLAGS)
	mv -f $@ $(QSSLIB_BIN_DIR) || exit 1

createMask3D:
	gcc $(DEBUG_FLAGS) $(QSSLIB_DIRS) \
	-I$(QSSLIB_INCLUDE) ./src/applications/createMask3D.c \
	-o $@ $(QSSLIB_LIBS) $(OMP_FLAGS)
	mv -f $@ $(QSSLIB_BIN_DIR) || exit 1

clean:
	cd src; make $@ || exit 1
	rm -f *.o
	rm -f *.tmp

cleanbin:
	cd $(QSSLIB_BIN_DIR); rm -f *

spotless: clean
	cd src; make $@ || exit 1
	cd examples; make $@ || exit 1
	if [ -d "include" ]; then rm -f -rf include/*; rmdir include; fi
	if [ -d "lib" ]; then rm -f -rf lib/*; rmdir lib; fi

