################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/pso/src/strategy/ibex_Optimizer_PSO.cpp \
../plugins/pso/src/strategy/ibex_TreeCellOpt.cpp 

OBJS += \
./plugins/pso/src/strategy/ibex_Optimizer_PSO.o \
./plugins/pso/src/strategy/ibex_TreeCellOpt.o 

CPP_DEPS += \
./plugins/pso/src/strategy/ibex_Optimizer_PSO.d \
./plugins/pso/src/strategy/ibex_TreeCellOpt.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/pso/src/strategy/%.o: ../plugins/pso/src/strategy/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


