################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Galaxies/create_sersic.c 

OBJS += \
./Galaxies/create_sersic.o 

C_DEPS += \
./Galaxies/create_sersic.d 


# Each subdirectory must supply rules for building sources it contributes
Galaxies/%.o: ../Galaxies/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"/Users/mpetkova/Code/NR/include" -I"/Users/mpetkova/Code/CosmoLib/include" -I"/Users/mpetkova/Code/SLsimLib/include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


