SRC_DIR = .

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config


# Set MAKE
    
all:  progs

progs:  drain createInitFromMask imbibe createMaskUBC 

drain : 
	gcc $(OPTIM_FLAGS) $(QSSLIB_DIRS) \
	-I$(QSSLIB_INCLUDE) drain.c \
	-o $@ $(QSSLIB_LIBS) $(OMP_FLAGS)
	mv -f $@ $(QSSLIB_BIN_DIR);

imbibe : 
	gcc $(OPTIM_FLAGS) $(QSSLIB_DIRS) \
	-I$(QSSLIB_INCLUDE) imbibe.c \
	-o $@ $(QSSLIB_LIBS) $(OMP_FLAGS)
	mv -f $@ $(QSSLIB_BIN_DIR);

createInitFromMask:
	gcc $(OPTIM_FLAGS) $(QSSLIB_DIRS) \
	-I$(QSSLIB_INCLUDE) createInitFromMask.c \
	$(QSSLIB_LIBS)  -o $@
	mv -f $@ $(QSSLIB_BIN_DIR);

createMaskUBC:
	gcc $(OPTIM_FLAGS) $(QSSLIB_DIRS) \
	-I$(QSSLIB_INCLUDE) createMaskUBC.c \
	$(QSSLIB_LIBS)  -o $@
	mv -f $@ $(QSSLIB_BIN_DIR);

clean:
	$(RM) *.o
