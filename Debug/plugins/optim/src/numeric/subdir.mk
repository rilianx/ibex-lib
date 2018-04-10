################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/optim/src/numeric/ibex_LinearizerCombo.cpp 

OBJS += \
./plugins/optim/src/numeric/ibex_LinearizerCombo.o 

CPP_DEPS += \
./plugins/optim/src/numeric/ibex_LinearizerCombo.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/optim/src/numeric/%.o: ../plugins/optim/src/numeric/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


