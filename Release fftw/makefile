RM := rm -rf

-include src/subdir.mk

# All Target
all: libQuaSAR.so

# Tool invocations
libQuaSAR.so: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -shared -o "libQuaSAR.so" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(LIBRARIES)$(CC_DEPS)$(C++_DEPS)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) libQuaSAR.so
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
