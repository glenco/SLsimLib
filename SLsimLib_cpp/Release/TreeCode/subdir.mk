################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../TreeCode/TreeNB.c \
../TreeCode/TreeNBDriver.c \
../TreeCode/TreeNBForce.c \
../TreeCode/rayshooterNB.c \
../TreeCode/readfiles.c 

OBJS += \
./TreeCode/TreeNB.o \
./TreeCode/TreeNBDriver.o \
./TreeCode/TreeNBForce.o \
./TreeCode/rayshooterNB.o \
./TreeCode/readfiles.o 

C_DEPS += \
./TreeCode/TreeNB.d \
./TreeCode/TreeNBDriver.d \
./TreeCode/TreeNBForce.d \
./TreeCode/rayshooterNB.d \
./TreeCode/readfiles.d 


# Each subdirectory must supply rules for building sources it contributes
TreeCode/%.o: ../TreeCode/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"/Users/mpetkova/Code/NR/include" -I"/Users/mpetkova/Code/CosmoLib/include" -I"/Users/mpetkova/Code/SLsimLib/include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


