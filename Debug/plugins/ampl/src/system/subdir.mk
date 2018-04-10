################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/ampl/src/system/ibex_AmplInterface.cpp 

OBJS += \
./plugins/ampl/src/system/ibex_AmplInterface.o 

CPP_DEPS += \
./plugins/ampl/src/system/ibex_AmplInterface.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/ampl/src/system/%.o: ../plugins/ampl/src/system/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


