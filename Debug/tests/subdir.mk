################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tests/ExFunction.cpp \
../tests/Instance.cpp \
../tests/Ponts30.cpp \
../tests/TestAgenda.cpp \
../tests/TestArith.cpp \
../tests/TestBitSet.cpp \
../tests/TestBoolInterval.cpp \
../tests/TestCtcExist.cpp \
../tests/TestCtcForAll.cpp \
../tests/TestCtcFwdBwd.cpp \
../tests/TestCtcHC4.cpp \
../tests/TestCtcInteger.cpp \
../tests/TestCtcNotIn.cpp \
../tests/TestDim.cpp \
../tests/TestDomain.cpp \
../tests/TestDoubleIndex.cpp \
../tests/TestEval.cpp \
../tests/TestExpr.cpp \
../tests/TestExpr2DAG.cpp \
../tests/TestExpr2Minibex.cpp \
../tests/TestExprCmp.cpp \
../tests/TestExprCopy.cpp \
../tests/TestExprDiff.cpp \
../tests/TestExprLinearity.cpp \
../tests/TestExprSimplify.cpp \
../tests/TestFritzJohn.cpp \
../tests/TestFunction.cpp \
../tests/TestGradient.cpp \
../tests/TestHC4Revise.cpp \
../tests/TestInHC4Revise.cpp \
../tests/TestInnerArith.cpp \
../tests/TestInterval.cpp \
../tests/TestIntervalMatrix.cpp \
../tests/TestIntervalVector.cpp \
../tests/TestLinear.cpp \
../tests/TestNewton.cpp \
../tests/TestNumConstraint.cpp \
../tests/TestParser.cpp \
../tests/TestPdcHansenFeasibility.cpp \
../tests/TestSeparator.cpp \
../tests/TestSet.cpp \
../tests/TestSolver.cpp \
../tests/TestString.cpp \
../tests/TestSymbolMap.cpp \
../tests/TestSystem.cpp \
../tests/TestSystemBox.cpp \
../tests/TestTimer.cpp \
../tests/TestVarSet.cpp \
../tests/utest.cpp \
../tests/utils.cpp 

OBJS += \
./tests/ExFunction.o \
./tests/Instance.o \
./tests/Ponts30.o \
./tests/TestAgenda.o \
./tests/TestArith.o \
./tests/TestBitSet.o \
./tests/TestBoolInterval.o \
./tests/TestCtcExist.o \
./tests/TestCtcForAll.o \
./tests/TestCtcFwdBwd.o \
./tests/TestCtcHC4.o \
./tests/TestCtcInteger.o \
./tests/TestCtcNotIn.o \
./tests/TestDim.o \
./tests/TestDomain.o \
./tests/TestDoubleIndex.o \
./tests/TestEval.o \
./tests/TestExpr.o \
./tests/TestExpr2DAG.o \
./tests/TestExpr2Minibex.o \
./tests/TestExprCmp.o \
./tests/TestExprCopy.o \
./tests/TestExprDiff.o \
./tests/TestExprLinearity.o \
./tests/TestExprSimplify.o \
./tests/TestFritzJohn.o \
./tests/TestFunction.o \
./tests/TestGradient.o \
./tests/TestHC4Revise.o \
./tests/TestInHC4Revise.o \
./tests/TestInnerArith.o \
./tests/TestInterval.o \
./tests/TestIntervalMatrix.o \
./tests/TestIntervalVector.o \
./tests/TestLinear.o \
./tests/TestNewton.o \
./tests/TestNumConstraint.o \
./tests/TestParser.o \
./tests/TestPdcHansenFeasibility.o \
./tests/TestSeparator.o \
./tests/TestSet.o \
./tests/TestSolver.o \
./tests/TestString.o \
./tests/TestSymbolMap.o \
./tests/TestSystem.o \
./tests/TestSystemBox.o \
./tests/TestTimer.o \
./tests/TestVarSet.o \
./tests/utest.o \
./tests/utils.o 

CPP_DEPS += \
./tests/ExFunction.d \
./tests/Instance.d \
./tests/Ponts30.d \
./tests/TestAgenda.d \
./tests/TestArith.d \
./tests/TestBitSet.d \
./tests/TestBoolInterval.d \
./tests/TestCtcExist.d \
./tests/TestCtcForAll.d \
./tests/TestCtcFwdBwd.d \
./tests/TestCtcHC4.d \
./tests/TestCtcInteger.d \
./tests/TestCtcNotIn.d \
./tests/TestDim.d \
./tests/TestDomain.d \
./tests/TestDoubleIndex.d \
./tests/TestEval.d \
./tests/TestExpr.d \
./tests/TestExpr2DAG.d \
./tests/TestExpr2Minibex.d \
./tests/TestExprCmp.d \
./tests/TestExprCopy.d \
./tests/TestExprDiff.d \
./tests/TestExprLinearity.d \
./tests/TestExprSimplify.d \
./tests/TestFritzJohn.d \
./tests/TestFunction.d \
./tests/TestGradient.d \
./tests/TestHC4Revise.d \
./tests/TestInHC4Revise.d \
./tests/TestInnerArith.d \
./tests/TestInterval.d \
./tests/TestIntervalMatrix.d \
./tests/TestIntervalVector.d \
./tests/TestLinear.d \
./tests/TestNewton.d \
./tests/TestNumConstraint.d \
./tests/TestParser.d \
./tests/TestPdcHansenFeasibility.d \
./tests/TestSeparator.d \
./tests/TestSet.d \
./tests/TestSolver.d \
./tests/TestString.d \
./tests/TestSymbolMap.d \
./tests/TestSystem.d \
./tests/TestSystemBox.d \
./tests/TestTimer.d \
./tests/TestVarSet.d \
./tests/utest.d \
./tests/utils.d 


# Each subdirectory must supply rules for building sources it contributes
tests/%.o: ../tests/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


