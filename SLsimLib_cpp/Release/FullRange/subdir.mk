################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../FullRange/implant_stars.c \
../FullRange/internal_rayshooter_nfw.c 

OBJS += \
./FullRange/implant_stars.o \
./FullRange/internal_rayshooter_nfw.o 

C_DEPS += \
./FullRange/implant_stars.d \
./FullRange/internal_rayshooter_nfw.d 


# Each subdirectory must supply rules for building sources it contributes
FullRange/%.o: ../FullRange/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"/Users/mpetkova/Code/NR/include" -I"/Users/mpetkova/Code/CosmoLib/include" -I"/Users/mpetkova/Code/SLsimLib/include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


