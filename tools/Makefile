#--gcc
MPI_CC         = mpicc
MPI_F77        = mpif90
MPI_F90        = mpif90
CC         = gcc
F77        = gfortran
F90        = gfortran
#--intel fortran
#MPI_CC         = mpiicc
#MPI_F77        = mpiifort
#MPI_F90        = mpiifort
#CC         = icc
#F77        = ifort
#F90        = ifort

DEBUG_FLAGS=-g -fbacktrace -gdwarf-2 -fbounds-check -finit-real=nan -finit-integer=nan -ffpe-trap=invalid,zero,overflow

#--============================== SUPER LU -=====================================#
SUPERLU_LIBS=-L/data1/ghosh/installations/SuperLU_DIST_4.1/lib -lsuperlu_dist_4.1
SUPERLU_INC=-I/data1/ghosh/installations/SuperLU_DIST_4.1/SRC
#============================ BLAS and LAPACK ===================================#
#--(A) use intel-MKL for blas+lapack (assuming env var has been set)
BLAS_LIBS=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#--(B) use gcc blas+lapack on the server
#BLAS_LIBS= -L/usr/lib64/ -llapack -lblas
#============================ METIS and PARMETIS ================================#
#--(A) use metis+parmetis, built with gcc-4.9.1 and OpenMPI-1.4.3 for 64bit integers and 64bit floats -- metis gives segfault
#METIS_LIB=/data1/ghosh/installations/parmetis-4.0.3/build/Linux-x86_64/libmetis/libmetis.a
#PARMETIS_LIB=/data1/ghosh/installations/parmetis-4.0.3/build/Linux-x86_64/libparmetis/libparmetis.a
#--(B) use HHPC-default metis+parmetis at /use/local/lib
METIS_LIB=-lmetis
PARMETIS_LIB=-lparmetis
#--------------------------------------------------------------------------------#
METISPARMETIS_LIBS=${PARMETIS_LIB} ${METIS_LIB}
#================================================================================#

all: yieldStress
all: yieldStressCMRL
all: tpConv
all: tpGrains
all: scaleMesh
all: calcSchmidFactorsHCP
all: genTexture
all: plotTexture
all: matchCrystallography
all: calcPatches
all: insertFieldIntoGrain
all: convNye
all: calcPlasticStrain
all: calcEqvPlasticStrainMacroscopic
all: calcNyeNorm
all: nodalRecoveryProg
all: regularizeField
all: regularizeFieldParallel
all: printGrainNeighbors
all: calcCorrelation
all: findFieldLocalMaxima
all: calculateFieldStatistics
all: findCriticalHardSoftPairs
all: convABQ
all: genInp
all: eulerToSchmidFactors
all: convertProcessorBasedOutputs
all: convertProcessorBasedOutputsBB
all: convTextureToCMRL
all: correctConventionCMRL

yieldStress: yieldStress.f90
	gfortran yieldStress.f90 -g -O2 -o yieldStress
yieldStressCMRL: yieldStressCMRL.f90
	gfortran -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow yieldStressCMRL.f90 -g -O2 -o yieldStressCMRL
tpConv: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 tp.f90
	gfortran ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 tp.f90 -g -O3 -o tpConv
plotElementResponse: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 plotElementResponse.f90
	gfortran ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 plotElementResponse.f90 -g -O3 -o plotElementResponse
tpGrains: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/crystal_Titanium.f90 tpGrains.f90
	gfortran ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/crystal_Titanium.f90 tpGrains.f90 ${DEBUG_FLAGS} -O0 -o tpGrains
scaleMesh: scaleMesh.90
	gfortran scaleMesh.f90  -fbacktrace -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow -O3  -o scaleMesh
calcPlasticStrain: calcPlasticStrain.f90
	gfortran calcPlasticStrain.f90 ${DEBUG_FLAGS} -O3 -o calcPlasticStrain
calcEqvPlasticStrainMacroscopic: calcEqvPlasticStrainMacroscopic.f90 ../../HomogenizedModelTi/functions/fortran/stressAnalysis.f90
	gfortran calcEqvPlasticStrainMacroscopic.f90 ../../HomogenizedModelTi/functions/fortran/stressAnalysis.f90 ${DEBUG_FLAGS} -O3 -o calcEqvPlasticStrainMacroscopic
genTexture: ../modules/strings.f90 ../modules/init_random_seed.f90 ./genTexture.f90
	gfortran ../modules/strings.f90 ../modules/init_random_seed.f90 ./genTexture.f90 -g -O3 -o ./genTexture
plotTexture: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../src/options_mod.f90 ../modules/crystal_Titanium.f90 ../modules/strings.f90 ./plotTexture.f90
	gfortran ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../src/options_mod.f90 ../modules/crystal_Titanium.f90 ../modules/strings.f90 ./plotTexture.f90 -g -O3 -o ./plotTexture
matchCrystallography: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/crystal_Titanium.f90 ../../HomogenizedModelTi/functions/fortran/calculateMisorientationAngle.f90 ../../HomogenizedModelTi/functions/fortran/findContainerInterval.f90 ../../HomogenizedModelTi/functions/fortran/readGrainNeighborListsCPFE.f90 ../../HomogenizedModelTi/functions/fortran/getRotationMatrix.f90 ../../HomogenizedModelTi/functions/fortran/init_random_seed.f90  ../../HomogenizedModelTi/functions/fortran/calculateBoxMisorientationHistogram.f90 matchCrystallography.f90
	gfortran ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/crystal_Titanium.f90 ../../HomogenizedModelTi/functions/fortran/calculateMisorientationAngle.f90 ../../HomogenizedModelTi/functions/fortran/findContainerInterval.f90 ../../HomogenizedModelTi/functions/fortran/readGrainNeighborListsCPFE.f90 ../../HomogenizedModelTi/functions/fortran/getRotationMatrix.f90 ../../HomogenizedModelTi/functions/fortran/init_random_seed.f90  ../../HomogenizedModelTi/functions/fortran/calculateBoxMisorientationHistogram.f90 matchCrystallography.f90 -g -O3 -o matchCrystallography
shuffleCrystallography: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/crystal_Titanium.f90 ../../HomogenizedModelTi/functions/fortran/sortNumbers.f90 ../../HomogenizedModelTi/functions/fortran/saveGrainTexture.f90 ../../HomogenizedModelTi/functions/fortran/init_random_seed.f90 shuffleCrystallography.f90
	gfortran  ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/crystal_Titanium.f90 ../../HomogenizedModelTi/functions/fortran/sortNumbers.f90 ../../HomogenizedModelTi/functions/fortran/saveGrainTexture.f90 ../../HomogenizedModelTi/functions/fortran/init_random_seed.f90 shuffleCrystallography.f90 -g -O3 -o shuffleCrystallography
calcSchmidFactorsHCP: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/crystal_Titanium.f90 ../modules/strings.f90 ./calcSchmidFactorsHCP.f90
	gfortran ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/crystal_Titanium.f90 ../modules/strings.f90 ./calcSchmidFactorsHCP.f90 -g -O0 -o ./calcSchmidFactorsHCP
calcPatches: ../modules/meshtools_mod.f90 ../modules/nodalRecovery_mod.f90 ../modules/grains_mod.f90 ./calcPatches.f90
	${F90} ../modules/meshtools_mod.f90 ../modules/nodalRecovery_mod.f90 ../modules/grains_mod.f90 ./calcPatches.f90 ${BLAS_LIBS} -o calcPatches
insertFieldIntoGrain:
	${F90} ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 insertUniformFieldIntoGrain.f90 -o insertUniformFieldIntoGrain
convNye:	../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/nodalRecovery_mod.f90 ../modules/strings.f90 ../src/options_mod.f90 ./convNye.f90
	${MPI_F90} ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/nodalRecovery_mod.f90 ../modules/strings.f90 ../src/options_mod.f90 \
			   convNye.f90 ${BLAS_LIBS} ${DEBUG_FLAGS} -O3 -o convNye
calcNyeNorm: calcNyeNorm.f90
	${F90} calcNyeNorm.f90 -O3 -o calcNyeNorm
calcTriaxiality: calcTriaxiality.f90
	${F90} calcTriaxiality.f90 ${DEBUG_FLAGS} -O2 -o calcTriaxiality
nodalRecoveryProg: ../modules/meshtools_mod.f90 ../modules/nodalRecovery_mod.f90 ../modules/grains_mod.f90 nodalRecoveryProg.f90
	${F90} ../modules/meshtools_mod.f90 ../modules/nodalRecovery_mod.f90 ../modules/grains_mod.f90 \
			   nodalRecoveryProg.f90 ${DEBUG_FLAGS} ${BLAS_LIBS} -O3 -o nodalRecoveryProg
regularizeField: ../modules/strings.f90 ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/regularization_mod.f90 regularizeField.f90
	${F90} ../modules/strings.f90 ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/regularization_mod.f90 \
			   regularizeField.f90 ${DEBUG_FLAGS} -O3 -o regularizeField
regularizeFieldParallel:  ../modules/strings.f90 ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/regularization_mod.f90 regularizeFieldParallel.f90
	${MPI_F90} ../modules/strings.f90 ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/regularization_mod.f90 \
			   regularizeFieldParallel.f90 -g -fbacktrace -gdwarf-2 -O3 -fbounds-check -finit-real=nan -finit-integer=nan -ffpe-trap=invalid,zero,overflow -o regularizeFieldParallel
printGrainNeighbors: ../modules/meshtools_mod.f90 ../modules/strings.f90 ../modules/grains_mod.f90 printGrainNeighbors.f90
	${F90} ../modules/meshtools_mod.f90 ../modules/strings.f90 ../modules/grains_mod.f90 printGrainNeighbors.f90 -o printGrainNeighbors
calcCorrelation: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 calcCorrelation.f90
	${F90} ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 calcCorrelation.f90 -g -O0 -fbacktrace -gdwarf-2 -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wuninitialized -o calcCorrelation
findFieldLocalMaxima: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 findFieldLocalMaxima.f90
	${F90} ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 findFieldLocalMaxima.f90 -g -O0 -gdwarf-2 -fbacktrace -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wuninitialized -o findFieldLocalMaxima
calculateFieldStatistics: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../../HomogenizedModelTi/functions/fortran/stressAnalysis.f90 calculateFieldStatistics.f90
	${F90} ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../../HomogenizedModelTi/functions/fortran/stressAnalysis.f90 calculateFieldStatistics.f90 -g -O3 -gdwarf-2 -fbacktrace -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wuninitialized -o calculateFieldStatistics
findCriticalHardSoftPairs: ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 findCriticalHardSoftPairs.f90
	${F90} ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 findCriticalHardSoftPairs.f90 -g -O0 -gdwarf-2 -fbacktrace -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wuninitialized -o findCriticalHardSoftPairs
convABQ: ../modules/meshtools_mod.f90 ../modules/init_random_seed.f90 convABQ.f90
	${F90} ../modules/meshtools_mod.f90 ../modules/init_random_seed.f90 convABQ.f90 -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow ${BLAS_LIBS} -g -O2 -o convABQ
genInp: genInp.f90 ../modules/init_random_seed.f90
	${F90} genInp.f90 ../modules/init_random_seed.f90 -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow -g -O0 -o genInp
eulerToSchmidFactors: ../modules/crystal_Titanium.f90 eulerToSchmidFactors.f90
	gfortran ../modules/crystal_Titanium.f90 ./eulerToSchmidFactors.f90 -g -fbacktrace -gdwarf-2 -fbounds-check -finit-real=nan -finit-integer=nan -ffpe-trap=invalid,zero,overflow -O0 -o eulerToSchmidFactors
interpolateTimeStep: interpolateTimeStep.f90
	gfortran interpolateTimeStep.f90 -g -fbacktrace -gdwarf-2 -O3 -fbounds-check -finit-real=nan -finit-integer=nan -ffpe-trap=invalid,zero,overflow -o interpolateTimeStep
convertProcessorBasedOutputsBB: convertProcessorBasedOutputsBB.f90
	gfortran convertProcessorBasedOutputsBB.f90 -g -fbacktrace -gdwarf-2 -fbounds-check -finit-real=nan -finit-integer=nan -ffpe-trap=invalid,zero,overflow -O3 -o convertProcessorBasedOutputsBB
convertProcessorBasedOutputs: convertProcessorBasedOutputs.f90
	gfortran convertProcessorBasedOutputs.f90 -g -fbacktrace -gdwarf-2 -fbounds-check -finit-real=nan -finit-integer=nan -ffpe-trap=invalid,zero,overflow -O3 -o convertProcessorBasedOutputs
convTextureToCMRL: ../modules/meshtools_mod.f90 ../modules/crystal_Titanium.f90 ../modules/grains_mod.f90 convTextureToCMRL.f90
	gfortran ../modules/meshtools_mod.f90 ../modules/grains_mod.f90 ../modules/crystal_Titanium.f90 convTextureToCMRL.f90 -g -fbacktrace -gdwarf-2 -fbounds-check -finit-real=nan -finit-integer=nan -ffpe-trap=invalid,zero,overflow -O3 -o convTextureToCMRL
correctConventionCMRL: correctConventionCMRL.f90
	gfortran  correctConventionCMRL.f90 ${DEBUG_FLAGS} -O0 -o correctConventionCMRL
