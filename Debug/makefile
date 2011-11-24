################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include TreeCode_link/List/subdir.mk
-include TreeCode_link/subdir.mk
-include TreeCode/subdir.mk
-include Kist/subdir.mk
-include ImageProcessing/subdir.mk
-include Galaxies/subdir.mk
-include FullRange/subdir.mk
-include Fitlens/subdir.mk
-include BLR/subdir.mk
-include AnalyticNSIE/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: libSLsimLib_cpp.a

# Tool invocations
libSLsimLib_cpp.a: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC Archiver'
	ar s -r  "libSLsimLib_cpp.a" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(C++_DEPS)$(OBJS)$(C_DEPS)$(CC_DEPS)$(ARCHIVES)$(CPP_DEPS)$(CXX_DEPS)$(C_UPPER_DEPS) libSLsimLib_cpp.a
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
