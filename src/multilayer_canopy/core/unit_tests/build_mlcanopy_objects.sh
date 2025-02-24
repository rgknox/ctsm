#!/bin/bash

# Build Multilayer Canopy F90 Objects

FC='gfortran'

F_OPTS="-shared -fPIC" # -O3"

#F_OPTS="-shared -fPIC -O0 -g -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -Wall"

MOD_FLAG="-J"

rm -f bld/*.o
rm -f bld/*.mod

# Build the new file with constants

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/shr_kind_mod.o share/shr_kind_mod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/shr_const_mod.o share/shr_const_mod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/shr_sys_mod.o share/shr_sys_mod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLCanopyVarPar.o ../src/MLCanopyVarPar.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLCanopyVarCon.o ../src/MLCanopyVarCon.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLCanopyVarCtl.o ../src/MLCanopyVarCtl.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/SimpleGetPsiHatMod.o share/SimpleGetPsiHatMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLCanopyFluxesType.o ../src/MLCanopyFluxesType.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLMathToolsMod.o ../src/MLMathToolsMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLWaterVaporMod.o ../src/MLWaterVaporMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MlLeafPhotosynthesisMod.o ../src/MLLeafPhotosynthesisMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLSolarRadiationMod.o ../src/MLSolarRadiationMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLCanopyNitrogenProfileMod.o ../src/MLCanopyNitrogenProfileMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLCanopyWaterMod.o ../src/MLCanopyWaterMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLLongwaveRadiationMod.o ../src/MLLongwaveRadiationMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLWaterVaporMod.o ../src/MLWaterVaporMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLLeafFluxesMod.o ../src/MLLeafFluxesMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLSoilFluxesMod.o ../src/MLSoilFluxesMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLCanopyTurbulenceMod.o ../src/MLCanopyTurbulenceMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLLeafHeatCapacityMod.o ../src/MLLeafHeatCapacityMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLPlantHydraulicsMod.o ../src/MLPlantHydraulicsMod.F90
#${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/



