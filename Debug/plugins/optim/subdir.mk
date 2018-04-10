################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/optim/benchmark_optim.cpp 

OBJS += \
./plugins/optim/benchmark_optim.o 

CPP_DEPS += \
./plugins/optim/benchmark_optim.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/optim/%.o: ../plugins/optim/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


