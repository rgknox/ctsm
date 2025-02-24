#!/bin/bash

# Build Multilayer Canopy F90 Objects

FC='gfortran'

F_OPTS="-shared -fPIC -O3"

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
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/MLCanopyFluxesType.o ../src/MLCanopyFluxesType.F90

