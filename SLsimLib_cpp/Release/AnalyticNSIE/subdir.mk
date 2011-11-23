################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../AnalyticNSIE/expanded_powerlaw.c \
../AnalyticNSIE/lens_expand.c \
../AnalyticNSIE/mark_points.c \
../AnalyticNSIE/nfw_lens.c \
../AnalyticNSIE/nsie.c \
../AnalyticNSIE/powerlaw.c \
../AnalyticNSIE/randomize_lens.c \
../AnalyticNSIE/readfiles_ana.c 

OBJS += \
./AnalyticNSIE/expanded_powerlaw.o \
./AnalyticNSIE/lens_expand.o \
./AnalyticNSIE/mark_points.o \
./AnalyticNSIE/nfw_lens.o \
./AnalyticNSIE/nsie.o \
./AnalyticNSIE/powerlaw.o \
./AnalyticNSIE/randomize_lens.o \
./AnalyticNSIE/readfiles_ana.o 

C_DEPS += \
./AnalyticNSIE/expanded_powerlaw.d \
./AnalyticNSIE/lens_expand.d \
./AnalyticNSIE/mark_points.d \
./AnalyticNSIE/nfw_lens.d \
./AnalyticNSIE/nsie.d \
./AnalyticNSIE/powerlaw.d \
./AnalyticNSIE/randomize_lens.d \
./AnalyticNSIE/readfiles_ana.d 


# Each subdirectory must supply rules for building sources it contributes
AnalyticNSIE/%.o: ../AnalyticNSIE/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"/Users/mpetkova/Code/NR/include" -I"/Users/mpetkova/Code/CosmoLib/include" -I"/Users/mpetkova/Code/SLsimLib/include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


