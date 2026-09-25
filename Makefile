CC=g++
CFLAGS=-c -g -Wall -O3 -pthread
# Default EOS table directory baked into the binary (overridable at runtime with $TWOFLUID_EOS_DIR)
CPPFLAGS=-DTWOFLUID_DEFAULT_EOS_DIR=\"$(CURDIR)/eos_tables\"
SRC_DIR=src
OBJ_DIR=build
HEAD_DIR=include

EXEC=TWO_FLUID_TOV

SOURCES=$(addprefix $(SRC_DIR)/, \
	main.cpp conversions.cpp eos_analytic.cpp eos_tabular.cpp eos_poly.cpp eos_SIDM.cpp twofluid_TOV.cpp TOV_family.cpp )

OBJECTS=$(SOURCES:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CC) $^ -o $@ -lgsl -lgslcblas -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXEC)
