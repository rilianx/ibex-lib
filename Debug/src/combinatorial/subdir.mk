################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/combinatorial/ibex_QInter.cpp 

OBJS += \
./src/combinatorial/ibex_QInter.o 

CPP_DEPS += \
./src/combinatorial/ibex_QInter.d 


# Each subdirectory must supply rules for building sources it contributes
src/combinatorial/%.o: ../src/combinatorial/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


