################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../TreeCode_link/List/List.cpp 

OBJS += \
./TreeCode_link/List/List.o 

CPP_DEPS += \
./TreeCode_link/List/List.d 


# Each subdirectory must supply rules for building sources it contributes
TreeCode_link/List/%.o: ../TreeCode_link/List/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/mpetkova/Code/NR/include" -I"/Users/mpetkova/Code/SLsimLib_cpp/include" -I"/Users/mpetkova/Code/CosmoLib/include" -I"/Users/mpetkova/Code/SimpleTree" -O3 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


