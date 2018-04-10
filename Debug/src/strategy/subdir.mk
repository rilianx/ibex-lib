################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/strategy/ibex_Cell.cpp \
../src/strategy/ibex_CellBuffer.cpp \
../src/strategy/ibex_CellList.cpp \
../src/strategy/ibex_CellStack.cpp \
../src/strategy/ibex_Paver.cpp \
../src/strategy/ibex_SetImage.cpp 

OBJS += \
./src/strategy/ibex_Cell.o \
./src/strategy/ibex_CellBuffer.o \
./src/strategy/ibex_CellList.o \
./src/strategy/ibex_CellStack.o \
./src/strategy/ibex_Paver.o \
./src/strategy/ibex_SetImage.o 

CPP_DEPS += \
./src/strategy/ibex_Cell.d \
./src/strategy/ibex_CellBuffer.d \
./src/strategy/ibex_CellList.d \
./src/strategy/ibex_CellStack.d \
./src/strategy/ibex_Paver.d \
./src/strategy/ibex_SetImage.d 


# Each subdirectory must supply rules for building sources it contributes
src/strategy/%.o: ../src/strategy/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


