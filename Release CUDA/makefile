OS_SUFFIX := $(subst Linux,linux,$(subst Darwin/x86_64,darwin,$(shell uname -s)/$(shell uname -m)))
PYTHON_LDFLAGS := $(shell python3.12-config --ldflags)

RM := rm -rf

-include src/subdir.mk

GENCODE = \
-gencode arch=compute_89,code=sm_89 

# OpenGL+GLUT OS-specific define
ifeq ($(shell uname -s),Darwin)
GLUT_LIBS := -Xlinker -framework -Xlinker GLUT -Xlinker -framework -Xlinker OpenGL
else
GLUT_LIBS := -lGL -lGLU -lglut 
endif


# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: libQuaSAR.so

# Tool invocations
libQuaSAR.so: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: NVCC Linker'
	/usr/local/cuda-12.6/bin/nvcc --cudart static -shared --relocatable-device-code=false $(GENCODE) -link -o  "libQuaSAR.so" $(OBJS) $(USER_OBJS) $(LIBS) $(PYTHON_LDFLAGS)

	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(LIBRARIES)$(CC_DEPS)$(C++_DEPS)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CU_DEPS)$(CPP_DEPS)$(C_DEPS) libQuaSAR.so
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
