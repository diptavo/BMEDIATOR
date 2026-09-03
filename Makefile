CXX      ?= g++
USE_OPENMP ?= 0
OPENMP_ROOT ?=

BASE_CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -I. -Wno-unused-variable -Wno-unused-function
BASE_LDFLAGS  = -lm

ifeq ($(USE_OPENMP),1)
  COMPILER_VERSION := $(shell $(CXX) --version 2>/dev/null)
  ifneq ($(findstring Apple clang,$(COMPILER_VERSION)),)
    ifeq ($(strip $(OPENMP_ROOT)),)
      OPENMP_ROOT := $(shell \
        if command -v brew >/dev/null 2>&1; then brew --prefix libomp 2>/dev/null; \
        elif test -d /opt/homebrew/opt/libomp; then printf '%s' /opt/homebrew/opt/libomp; \
        elif test -d /usr/local/opt/libomp; then printf '%s' /usr/local/opt/libomp; \
        fi)
    endif
    ifeq ($(strip $(OPENMP_ROOT)),)
      $(error USE_OPENMP=1 with Apple Clang requires libomp. Install Homebrew and run 'brew install libomp', or set OPENMP_ROOT=/path/to/libomp)
    endif
    BASE_CXXFLAGS += -Xpreprocessor -fopenmp -I$(OPENMP_ROOT)/include
    BASE_LDFLAGS  += -L$(OPENMP_ROOT)/lib -Wl,-rpath,$(OPENMP_ROOT)/lib -lomp
  else
    BASE_CXXFLAGS += -fopenmp
    BASE_LDFLAGS  += -fopenmp
  endif
endif

CXXFLAGS ?= $(BASE_CXXFLAGS)
LDFLAGS  ?= $(BASE_LDFLAGS)

BIN       = bmediator
BUILD_DIR = build
HDRS      = bmediator.h plink_ld.h gsmr_qc.h regional_ld.h
SRCS      = main.cpp io.cpp cavi.cpp factor_model.cpp plink_ld.cpp pipeline.cpp gsmr_qc.cpp regional_ld.cpp
OBJS      = $(BUILD_DIR)/main.o \
            $(BUILD_DIR)/io.o \
            $(BUILD_DIR)/cavi.o \
            $(BUILD_DIR)/factor_model.o \
            $(BUILD_DIR)/plink_ld.o \
            $(BUILD_DIR)/pipeline.o \
            $(BUILD_DIR)/gsmr_qc.o \
            $(BUILD_DIR)/regional_ld.o

.PHONY: all clean test test-regional-stress help

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

test-regional-stress: $(BIN)
	python3 sim/run_regional_ld_stress.py --binary ./$(BIN) --outdir build/regional_ld_stress --replicates 3

help:
	@echo "Targets:"
	@echo "  make              Build bmediator"
	@echo "  make test         Run the smoke test"
	@echo "  make test-regional-stress  Run the stochastic genotype/LD smoke stress"
	@echo "  make clean        Remove build artifacts"
	@echo ""
	@echo "Variables:"
	@echo "  USE_OPENMP=1      Enable OpenMP; Apple Clang requires Homebrew libomp"
	@echo "  OPENMP_ROOT=<dir> Override the detected libomp installation prefix"
	@echo "  CXX=<compiler>    Override C++ compiler"
