###########################################################

## USER SPECIFIC DIRECTORIES ##

# CUDA directory:
#CUDA_ROOT_DIR=/usr/local/cuda
CUDA_ROOT_DIR=/usr
#CUDA_ROOT_DIR=/cm/shared/apps/cuda11.2/toolkit/11.2.2
##########################################################

## CC COMPILER OPTIONS ##

# CC compiler options:
#CC := g++ -std=c++17 
CC := mpicxx -std=c++17 -g 
#CC_FLAGS := -L${HOME}/.local/lib -L/cm/shared/apps/geos/3.9.1/lib -lgeos -pthread
CC_FLAGS := -L/usr/lib/x86_64-linux-gnu -lgeos -lpthread
CC_LIBS :=

##########################################################

## NVCC COMPILER OPTIONS ##

# NVCC compiler options:
NVCC := nvcc -std=c++17 #-g -G  
NVCC_FLAGS := -w -m64 -gencode arch=compute_89,code=sm_89 #-gencode arch=compute_70,code=sm_70
NVCC_FLAGS_DM := $(NVCC_FLAGS) -rdc=true 
NVCC_LIBS :=

# CUDA library directory:
CUDA_LIB_DIR := -L$(CUDA_ROOT_DIR)/local/cuda-11.8/targets/x86_64-linux/lib
# CUDA include directory:
CUDA_INC_DIR := -I$(CUDA_ROOT_DIR)/local/cuda-11.8/targets/x86_64-linux/include
# CUDA linking libraries:
CUDA_LINK_LIBS := -lcudart -lcudadevrt 

##########################################################

## Project file structure ##
# Root directory:
MAKE_DIR = $(PWD)

# Source file directory:
SRC_DIR = $(MAKE_DIR)/src

# Main file directory:
MAIN_DIR = $(MAKE_DIR)/tests

# Object file directory:
OBJ_DIR = $(MAKE_DIR)/bin

# Target diretory:
APP_DIR  = $(MAKE_DIR)

# Header files:
INCLUDE := -isystem -I/usr/include -Iinclude $(CUDA_INC_DIR) 

##########################################################

## Make variables ##

# Target executable name:
TARGET = prog

# Source files:
CU_DYNAMIC_SRCS = $(SRC_DIR)/cuda_quadtree.cu
CU_SRCS := $(wildcard $(SRC_DIR)/*.cu)
CU_SRCS := $(filter-out $(CU_DYNAMIC_SRCS), $(CU_SRCS))
CPP_SRCS = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(MAIN_DIR)/*.cpp)

# Object files:
CU_DYNAMIC_OBJS = $(OBJ_DIR)/cuda_quadtree.o
CU_OBJS = $(addprefix $(OBJ_DIR)/,$(notdir $(CU_SRCS:.cu=.o)))
CPP_OBJS = $(addprefix $(OBJ_DIR)/,$(notdir $(CPP_SRCS:.cpp=.o)))

OBJS =  $(CPP_OBJS) $(CU_OBJS)  $(CU_DYNAMIC_OBJS) $(OBJ_DIR)/dlink.o

##########################################################

## Compile ##

# Link c++ and CUDA compiled object files to target executable:
all: build $(APP_DIR)/$(TARGET)

$(APP_DIR)/$(TARGET) : $(OBJS)
	$(CC) $(INCLUDE) $(CUDA_INC_DIR) $(OBJS) -o $@ $(CC_FLAGS) $(CUDA_LIB_DIR) $(CUDA_LINK_LIBS)

# dlink for CUDA Dynamic Parallelism
$(OBJ_DIR)/dlink.o : $(CU_OBJS)
	$(NVCC) $(NVCC_FLAGS) -dlink $(CU_DYNAMIC_OBJS) $(CU_OBJS) $< -o $@

$(OBJ_DIR)/cuda_quadtree.o : $(SRC_DIR)/cuda_quadtree.cu
	$(NVCC) $(NVCC_FLAGS_DM) $(INCLUDE) -c $< -o $@ $(NVCC_LIBS)

# Compile CUDA source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cu
	$(NVCC) $(NVCC_FLAGS) $(INCLUDE) -c $< -o $@ $(NVCC_LIBS)

# Compile C++ source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CC) $(CC_FLAGS) $(INCLUDE) -c $< -o $@

$(OBJ_DIR)/%.o : $(MAIN_DIR)/%.cpp
	$(CC) $(CC_FLAGS) $(INCLUDE) -c $< -o $@

build:
	@mkdir -p $(OBJ_DIR)

# Clean objects in object directory.
clean:
	$(RM) bin/* *.o $(APP_DIR)/$(TARGET)
