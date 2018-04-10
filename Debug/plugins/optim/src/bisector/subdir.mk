################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/optim/src/bisector/ibex_LSmear.cpp 

OBJS += \
./plugins/optim/src/bisector/ibex_LSmear.o 

CPP_DEPS += \
./plugins/optim/src/bisector/ibex_LSmear.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/optim/src/bisector/%.o: ../plugins/optim/src/bisector/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


