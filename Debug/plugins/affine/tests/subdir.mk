################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/affine/tests/TestAffineArith.cpp \
../plugins/affine/tests/TestAffineEval.cpp 

OBJS += \
./plugins/affine/tests/TestAffineArith.o \
./plugins/affine/tests/TestAffineEval.o 

CPP_DEPS += \
./plugins/affine/tests/TestAffineArith.d \
./plugins/affine/tests/TestAffineEval.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/affine/tests/%.o: ../plugins/affine/tests/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


