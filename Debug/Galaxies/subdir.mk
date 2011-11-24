################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Galaxies/create_sersic.cpp 

OBJS += \
./Galaxies/create_sersic.o 

CPP_DEPS += \
./Galaxies/create_sersic.d 


# Each subdirectory must supply rules for building sources it contributes
Galaxies/%.o: ../Galaxies/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/mpetkova/Code/SimpleTree" -I"/Users/mpetkova/Code/SLsimLib_cpp/include" -I"/Users/mpetkova/Code/NR/include" -I"/Users/mpetkova/Code/CosmoLib_cpp/include" -O3 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


