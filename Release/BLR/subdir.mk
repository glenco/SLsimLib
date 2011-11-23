################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../BLR/blr_surface_brightness.c \
../BLR/blr_surface_brightness2.0.c 

OBJS += \
./BLR/blr_surface_brightness.o \
./BLR/blr_surface_brightness2.0.o 

C_DEPS += \
./BLR/blr_surface_brightness.d \
./BLR/blr_surface_brightness2.0.d 


# Each subdirectory must supply rules for building sources it contributes
BLR/%.o: ../BLR/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"/Users/mpetkova/Code/NR/include" -I"/Users/mpetkova/Code/CosmoLib/include" -I"/Users/mpetkova/Code/SLsimLib/include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


