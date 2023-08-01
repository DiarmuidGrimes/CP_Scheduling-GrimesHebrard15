
MAINDIR = $(shell pwd)/../..


OPTFLAGS = -O3 

COMPILFLAGS = -Wall \
-D_UNIX \


CCC = g++ $(OPTFLAGS) $(COMPILFLAGS)

BIN=./bin
SRC=./lib/src
MOD=./models
OBJ=./lib/obj
INC=./include

OBJDIR = ./obj
SRCDIR = ./src


PMODELS = $(wildcard $(MOD)/src/*.cpp)
MODAUX = $(PMODELS:.cpp=)
MODELS = $(patsubst $(MOD)/src/%, %, $(MODAUX))
BINS = $(patsubst %, $(BIN)/%, $(MODELS))
EXE = $(BINS:=;)

PINCSRC = $(wildcard $(INC)/*.h)
PLIBSRC = $(wildcard $(SRC)/*.cpp)
PLIBAUX = $(PLIBSRC:.cpp=.o)
PLIBOBJ = $(patsubst $(SRC)/%, $(OBJ)/%, $(PLIBAUX))


CFLAGS = -I$(INC) -I$(SRCDIR)

LFLAGS = -L$(OBJ) 


scheduler: lib $(BIN)/scheduler

clean : 
	rm -rf $(OBJ)/*.o $(OBJ)/*.a $(SRC)/*~ $(MOD)/obj/*.o $(MOD)/src/*~ $(MOD)/src/*/*~ $(INC)/*~ *~ $(BIN)/* 



# The library
lib: $(PLIBOBJ)


$(OBJ)/%.o: $(SRC)/%.cpp $(INC)/%.h 
	$(CCC) $(CFLAGS) -c   $< -o $@ 

# The examples
$(BIN)/%: $(MOD)/obj/%.o $(PLIBOBJ) $(BIN) 
	$(CCC) $(CFLAGS) $< $(PLIBOBJ) -lm -o $@ $(LFLAGS)

$(MOD)/obj/%.o: $(MOD)/src/%.cpp $(MOD)/obj
	$(CCC) -c $(CFLAGS) -o $@ $<

