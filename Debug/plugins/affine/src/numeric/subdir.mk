################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/affine/src/numeric/ibex_LinearizerAffine2.cpp 

OBJS += \
./plugins/affine/src/numeric/ibex_LinearizerAffine2.o 

CPP_DEPS += \
./plugins/affine/src/numeric/ibex_LinearizerAffine2.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/affine/src/numeric/%.o: ../plugins/affine/src/numeric/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


