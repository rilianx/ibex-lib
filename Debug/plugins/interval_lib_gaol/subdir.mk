################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/interval_lib_gaol/wrapper.cpp 

OBJS += \
./plugins/interval_lib_gaol/wrapper.o 

CPP_DEPS += \
./plugins/interval_lib_gaol/wrapper.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/interval_lib_gaol/%.o: ../plugins/interval_lib_gaol/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


