################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../tests/3rd/fmemopen/fmemopen.c 

OBJS += \
./tests/3rd/fmemopen/fmemopen.o 

C_DEPS += \
./tests/3rd/fmemopen/fmemopen.d 


# Each subdirectory must supply rules for building sources it contributes
tests/3rd/fmemopen/%.o: ../tests/3rd/fmemopen/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


