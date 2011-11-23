################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../TreeCode_link/KistDriver.cpp \
../TreeCode_link/Tree.cpp \
../TreeCode_link/TreeDriver.cpp \
../TreeCode_link/change_redshits.cpp \
../TreeCode_link/curve_routines.cpp \
../TreeCode_link/divide_images.cpp \
../TreeCode_link/double_sort.cpp \
../TreeCode_link/find_crit.cpp \
../TreeCode_link/grid_maintenance.cpp \
../TreeCode_link/image_finder.cpp \
../TreeCode_link/image_finder_kist.cpp \
../TreeCode_link/map_images.cpp \
../TreeCode_link/tree_maintenance.cpp \
../TreeCode_link/utilities.cpp 

C_SRCS += \
../TreeCode_link/grid_initialization.c 

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
./TreeCode_link/grid_maintenance.o \
./TreeCode_link/image_finder.o \
./TreeCode_link/image_finder_kist.o \
./TreeCode_link/map_images.o \
./TreeCode_link/tree_maintenance.o \
./TreeCode_link/utilities.o 

C_DEPS += \
./TreeCode_link/grid_initialization.d 

CPP_DEPS += \
./TreeCode_link/KistDriver.d \
./TreeCode_link/Tree.d \
./TreeCode_link/TreeDriver.d \
./TreeCode_link/change_redshits.d \
./TreeCode_link/curve_routines.d \
./TreeCode_link/divide_images.d \
./TreeCode_link/double_sort.d \
./TreeCode_link/find_crit.d \
./TreeCode_link/grid_maintenance.d \
./TreeCode_link/image_finder.d \
./TreeCode_link/image_finder_kist.d \
./TreeCode_link/map_images.d \
./TreeCode_link/tree_maintenance.d \
./TreeCode_link/utilities.d 


# Each subdirectory must supply rules for building sources it contributes
TreeCode_link/%.o: ../TreeCode_link/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/mpetkova/Code/NR/include" -I"/Users/mpetkova/Code/SLsimLib_cpp/include" -I"/Users/mpetkova/Code/CosmoLib/include" -I"/Users/mpetkova/Code/SimpleTree" -O3 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

TreeCode_link/%.o: ../TreeCode_link/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


