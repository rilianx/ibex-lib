################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/continuation/src/ibex_Cont.cpp \
../plugins/continuation/src/ibex_ContCell.cpp \
../plugins/continuation/src/ibex_ParametricProof.cpp 

OBJS += \
./plugins/continuation/src/ibex_Cont.o \
./plugins/continuation/src/ibex_ContCell.o \
./plugins/continuation/src/ibex_ParametricProof.o 

CPP_DEPS += \
./plugins/continuation/src/ibex_Cont.d \
./plugins/continuation/src/ibex_ContCell.d \
./plugins/continuation/src/ibex_ParametricProof.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/continuation/src/%.o: ../plugins/continuation/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


