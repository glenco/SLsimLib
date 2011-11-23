################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include TreeCode_link/subdir.mk
-include TreeCode_link/List/subdir.mk
-include TreeCode/subdir.mk
-include Kist/subdir.mk
-include ImageProcessing/subdir.mk
-include Galaxies/subdir.mk
-include FullRange/subdir.mk
-include Fitlens/subdir.mk
-include BLR/subdir.mk
-include AnalyticNSIE/subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: libSLsimLib.a

# Tool invocations
libSLsimLib.a: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC Archiver'
	ar -r "libSLsimLib.a" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(ARCHIVES) libSLsimLib.a
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
