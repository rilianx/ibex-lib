# Affine arithmetic, vendored from github.com/ibex-team/ibex-affine (see UPSTREAM).
# The sources are upstream's, untouched; this file is the only local addition.

target_sources (ibex PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine.h
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine2_fAF2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine2_fAF2.h
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine3_fAFFullI.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine3_fAFFullI.h
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineDomain.h
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineMain.h
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineVar.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineVar.h
  ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineVector.h
  ${CMAKE_CURRENT_SOURCE_DIR}/function/ibex_AffineEval.h
  ${CMAKE_CURRENT_SOURCE_DIR}/contractor/ibex_CtcAffineFwdBwd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/contractor/ibex_CtcAffineFwdBwd.h
  ${CMAKE_CURRENT_SOURCE_DIR}/numeric/ibex_LinearizerAffine2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/numeric/ibex_LinearizerAffine2.h
)

target_include_directories (ibex PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/arithmetic>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/function>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/contractor>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/numeric>")
