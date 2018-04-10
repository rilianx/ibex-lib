################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/ampl/tests/TestAmpl.cpp 

OBJS += \
./plugins/ampl/tests/TestAmpl.o 

CPP_DEPS += \
./plugins/ampl/tests/TestAmpl.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/ampl/tests/%.o: ../plugins/ampl/tests/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


