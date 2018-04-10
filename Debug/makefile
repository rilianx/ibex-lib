################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include tests/3rd/fmemopen/subdir.mk
-include tests/subdir.mk
-include src/tools/subdir.mk
-include src/system/subdir.mk
-include src/symbolic/subdir.mk
-include src/strategy/subdir.mk
-include src/set/subdir.mk
-include src/predicate/subdir.mk
-include src/parser/subdir.mk
-include src/numeric/subdir.mk
-include src/function/subdir.mk
-include src/contractor/subdir.mk
-include src/combinatorial/subdir.mk
-include src/bisector/subdir.mk
-include src/arithmetic/subdir.mk
-include plugins/solver/src/subdir.mk
-include plugins/solver/main/subdir.mk
-include plugins/pso/src/strategy/subdir.mk
-include plugins/pso/src/pso/subdir.mk
-include plugins/pso/main/subdir.mk
-include plugins/param-estim/tests/subdir.mk
-include plugins/param-estim/src/cliquer/subdir.mk
-include plugins/param-estim/src/subdir.mk
-include plugins/param-estim/examples/subdir.mk
-include plugins/optim/tests/subdir.mk
-include plugins/optim/src/strategy/subdir.mk
-include plugins/optim/src/numeric/subdir.mk
-include plugins/optim/src/loup/subdir.mk
-include plugins/optim/src/contractor/subdir.mk
-include plugins/optim/src/bisector/subdir.mk
-include plugins/optim/main/subdir.mk
-include plugins/optim/subdir.mk
-include plugins/lp_lib_soplex/subdir.mk
-include plugins/lp_lib_none/subdir.mk
-include plugins/lp_lib_cplex/subdir.mk
-include plugins/lp_lib_clp/subdir.mk
-include plugins/interval_lib_gaol/subdir.mk
-include plugins/interval_lib_filib/subdir.mk
-include plugins/interval_lib_direct/subdir.mk
-include plugins/interval_lib_bias/subdir.mk
-include plugins/ensta-robotics/tests/subdir.mk
-include plugins/ensta-robotics/src/subdir.mk
-include plugins/continuation/tests/subdir.mk
-include plugins/continuation/src/subdir.mk
-include plugins/continuation/examples/subdir.mk
-include plugins/ampl/tests/subdir.mk
-include plugins/ampl/src/system/subdir.mk
-include plugins/ampl/examples/subdir.mk
-include plugins/affine/tests/subdir.mk
-include plugins/affine/src/numeric/subdir.mk
-include plugins/affine/src/function/subdir.mk
-include plugins/affine/src/arithmetic/subdir.mk
-include plugins/affine/examples/subdir.mk
-include plugins/affine-extended/src/arithmetic/subdir.mk
-include plugins/affine-extended/examples/subdir.mk
-include examples/slam/subdir.mk
-include examples/subdir.mk
-include __build__/src/tools/subdir.mk
-include __build__/src/system/subdir.mk
-include __build__/src/symbolic/subdir.mk
-include __build__/src/strategy/subdir.mk
-include __build__/src/set/subdir.mk
-include __build__/src/predicate/subdir.mk
-include __build__/src/parser/subdir.mk
-include __build__/src/numeric/subdir.mk
-include __build__/src/function/subdir.mk
-include __build__/src/contractor/subdir.mk
-include __build__/src/combinatorial/subdir.mk
-include __build__/src/bisector/subdir.mk
-include __build__/src/arithmetic/subdir.mk
-include __build__/plugins/solver/src/subdir.mk
-include __build__/plugins/solver/main/subdir.mk
-include __build__/plugins/pso/src/strategy/subdir.mk
-include __build__/plugins/pso/src/pso/subdir.mk
-include __build__/plugins/pso/main/subdir.mk
-include __build__/plugins/optim/src/strategy/subdir.mk
-include __build__/plugins/optim/src/numeric/subdir.mk
-include __build__/plugins/optim/src/loup/subdir.mk
-include __build__/plugins/optim/src/contractor/subdir.mk
-include __build__/plugins/optim/src/bisector/subdir.mk
-include __build__/plugins/optim/main/subdir.mk
-include __build__/plugins/ampl/src/system/subdir.mk
-include __build__/plugins/affine/src/numeric/subdir.mk
-include __build__/plugins/affine/src/function/subdir.mk
-include __build__/plugins/affine/src/arithmetic/subdir.mk
-include __build__/3rd/soplex-1.7.1/src/subdir.mk
-include __build__/3rd/soplex-1.7.1/obj/O.linux.x86_64.gnu.opt/lib/subdir.mk
-include __build__/3rd/soplex-1.7.1/obj/O.linux.x86_64.gnu.opt/bin/subdir.mk
-include __build__/3rd/soplex-1.7.1/doc/subdir.mk
-include __build__/3rd/mathlib-2.1.0/src/subdir.mk
-include __build__/3rd/gaol-4.2.0/gdtoa/test/subdir.mk
-include __build__/3rd/gaol-4.2.0/gdtoa/subdir.mk
-include __build__/3rd/gaol-4.2.0/gaol/subdir.mk
-include __build__/3rd/gaol-4.2.0/examples/subdir.mk
-include __build__/3rd/gaol-4.2.0/check/subdir.mk
-include __build__/3rd/amplsolvers/sys.x86_64.Linux/subdir.mk
-include __build__/3rd/amplsolvers/examples/subdir.mk
-include __build__/3rd/amplsolvers/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: ibex-pso

# Tool invocations
ibex-pso: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++  -o "ibex-pso" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) ibex-pso
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
