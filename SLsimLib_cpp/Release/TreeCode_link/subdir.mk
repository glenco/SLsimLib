################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../TreeCode_link/KistDriver.c \
../TreeCode_link/Tree.c \
../TreeCode_link/TreeDriver.c \
../TreeCode_link/change_redshits.c \
../TreeCode_link/curve_routines.c \
../TreeCode_link/divide_images.c \
../TreeCode_link/double_sort.c \
../TreeCode_link/find_crit.c \
../TreeCode_link/grid_initialization.c \
../TreeCode_link/image_finder.c \
../TreeCode_link/image_finder_kist.c \
../TreeCode_link/map_images.c \
../TreeCode_link/peak_refinement.c \
../TreeCode_link/utilities.c 

OBJS += \
./TreeCode_link/KistDriver.o \
./TreeCode_link/Tree.o \
./TreeCode_link/TreeDriver.o \
./TreeCode_link/change_redshits.o \
./TreeCode_link/curve_routines.o \
./TreeCode_link/divide_images.o \
./TreeCode_link/double_sort.o \
./TreeCode_link/find_crit.o \
./TreeCode_link/grid_initialization.o \
./TreeCode_link/image_finder.o \
./TreeCode_link/image_finder_kist.o \
./TreeCode_link/map_images.o \
./TreeCode_link/peak_refinement.o \
./TreeCode_link/utilities.o 

C_DEPS += \
./TreeCode_link/KistDriver.d \
./TreeCode_link/Tree.d \
./TreeCode_link/TreeDriver.d \
./TreeCode_link/change_redshits.d \
./TreeCode_link/curve_routines.d \
./TreeCode_link/divide_images.d \
./TreeCode_link/double_sort.d \
./TreeCode_link/find_crit.d \
./TreeCode_link/grid_initialization.d \
./TreeCode_link/image_finder.d \
./TreeCode_link/image_finder_kist.d \
./TreeCode_link/map_images.d \
./TreeCode_link/peak_refinement.d \
./TreeCode_link/utilities.d 


# Each subdirectory must supply rules for building sources it contributes
TreeCode_link/%.o: ../TreeCode_link/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"/Users/mpetkova/Code/NR/include" -I"/Users/mpetkova/Code/CosmoLib/include" -I"/Users/mpetkova/Code/SLsimLib/include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


