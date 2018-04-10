################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/pso/main/ibexpso.cpp 

OBJS += \
./plugins/pso/main/ibexpso.o 

CPP_DEPS += \
./plugins/pso/main/ibexpso.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/pso/main/%.o: ../plugins/pso/main/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


