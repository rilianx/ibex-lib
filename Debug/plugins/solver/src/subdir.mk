################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../plugins/solver/src/ibex_DefaultSolver.cpp \
../plugins/solver/src/ibex_Manifold.cpp \
../plugins/solver/src/ibex_Solver.cpp \
../plugins/solver/src/ibex_SolverOutputBox.cpp 

OBJS += \
./plugins/solver/src/ibex_DefaultSolver.o \
./plugins/solver/src/ibex_Manifold.o \
./plugins/solver/src/ibex_Solver.o \
./plugins/solver/src/ibex_SolverOutputBox.o 

CPP_DEPS += \
./plugins/solver/src/ibex_DefaultSolver.d \
./plugins/solver/src/ibex_Manifold.d \
./plugins/solver/src/ibex_Solver.d \
./plugins/solver/src/ibex_SolverOutputBox.d 


# Each subdirectory must supply rules for building sources it contributes
plugins/solver/src/%.o: ../plugins/solver/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


