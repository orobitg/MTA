PROJECT = mta

# Environment 
SRC = ./src
OBJ = ./obj
BIN = ./bin

# Better do not change this

MPICC = gcc
#MPICC = mpicc
#MPICC = /Users/orobitg/Applications/mpich2/bin/mpicc
#MPICC=/home/oro/Aplicacions/mpich2/bin/mpicc
GCCFLAGS := -I$(INC) -I$(MPIEXT)

# Add any options you like
GCCFLAGS += -g
GCCFLAGS += -O3
#MPIFLAGS=-DMPI_FLAG
MPIFLAGS=-DSEQ_FLAG
GCCFLAGS += ${MPIFLAGS}

LIBS = -lm -lpthread

MKDIR=mkdir
CP=cp

OBJS = 	\
	utils.o \
	parameters_utils.o \
	sequence_utils.o \
	distance_matrix_utils.o \
	guide_tree_utils.o \
	alignment_utils.o \
	evaluation_utils.o \
	main.o
	  
	
vpath %.c $(SRC)
vpath %.o   $(OBJ)
vpath %     $(BIN)

.PHONY: all dir clean distclean

all: dir build

build: $(PROJECT)

%.o: $(SRC)/%.c
	$(MPICC) $(GCCFLAGS) -c $(SRC)/$(<F) -o $(OBJ)/$(@F)

%: %.o
	$(MPICC) $(LDFLAGS) $(addprefix $(OBJ)/,$(filter %.o,$(^F))) -o $(BIN)/$(@F)

dir:
	mkdir -p $(BIN); mkdir -p $(OBJ);
	mkdir -p $(BIN)/plugins

clean:
	rm -f $(OBJ)/*.o

distclean:
	rm -rf $(BIN)
	rm -rf $(OBJ)

$(PROJECT): $(OBJS)
	$(MPICC) $(LDFLAGS) $(addprefix $(OBJ)/,$(OBJS)) $(LIBS) -o $(BIN)/$(PROJECT)

