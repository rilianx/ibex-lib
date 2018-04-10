################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/lp_lib_none/ibex_LPWrapper.cpp 

OBJS += \
./plugins/lp_lib_none/ibex_LPWrapper.o 

CPP_DEPS += \
./plugins/lp_lib_none/ibex_LPWrapper.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/lp_lib_none/%.o: ../plugins/lp_lib_none/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


