CC=g++
CFLAGS=-c -g -Ofast -pthread -I ./
SRC_DIR=src
OBJ_DIR=build
HEAD_DIR=include

EXEC=TWO_FLUID_TOV

SOURCES=$(addprefix $(SRC_DIR)/, \
	main.cpp conversions.cpp eos.cpp eos_tabular.cpp eos_poly.cpp twofluid_TOV.cpp TOV_family.cpp )

OBJECTS=$(SOURCES:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

all: $(EXEC)

$(EXEC): $(OBJECTS) 
	$(CC) $^ -o $@ -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@ -lm

clean: 
	rm $(OBJECTS) $(EXEC)