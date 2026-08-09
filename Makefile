CXX      ?= g++
USE_OPENMP ?= 0

BASE_CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -I. -Wno-unused-variable -Wno-unused-function
BASE_LDFLAGS  = -lm

ifeq ($(USE_OPENMP),1)
  BASE_CXXFLAGS += -fopenmp
  BASE_LDFLAGS  += -fopenmp
endif

CXXFLAGS ?= $(BASE_CXXFLAGS)
LDFLAGS  ?= $(BASE_LDFLAGS)

BIN       = bmediator
BUILD_DIR = build
HDRS      = bmediator.h plink_ld.h gsmr_qc.h
SRCS      = main.cpp io.cpp cavi.cpp plink_ld.cpp pipeline.cpp gsmr_qc.cpp
OBJS      = $(BUILD_DIR)/main.o \
            $(BUILD_DIR)/io.o \
            $(BUILD_DIR)/cavi.o \
            $(BUILD_DIR)/plink_ld.o \
            $(BUILD_DIR)/pipeline.o \
            $(BUILD_DIR)/gsmr_qc.o

.PHONY: all clean test help

all: $(BIN)

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(BIN) $(OBJS) $(LDFLAGS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: %.cpp $(HDRS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(BIN)

test: $(BIN)
	bash test/run_test.sh

help:
	@echo "Targets:"
	@echo "  make              Build bmediator"
	@echo "  make test         Run the smoke test"
	@echo "  make clean        Remove build artifacts"
	@echo ""
	@echo "Variables:"
	@echo "  USE_OPENMP=1      Enable OpenMP flags for compilers that support -fopenmp"
	@echo "  CXX=<compiler>    Override C++ compiler"
