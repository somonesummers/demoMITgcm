#!/usr/bin/env bash
# Script for building MITgcmuv
# Paul Summers Sept 2025
# invoked as:
# bash makeBuild.sh [root dir] [-mpi]
# Root directory is assumed to be ../../..
# -mpi builds with for mpi running. This requires specific packages to be loaded, currently this is set up for GaTech Pace
# MPI_HOME will need to be adjusted for local configuration. 
set -e

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: assuming root location of ../../.."
  ROOT="../../.."
elif [ "$1" == "-mpi" ]; then
  ROOT="../../.."
  OPT2=$1 
else
  ROOT=$1
  OPT2=$2
fi


unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     MACHINE="Linux";;
    Darwin*)    MACHINE="Mac";;
esac
echo "Idenitfied machine as ${MACHINE}"

cd build
FILE=Makefile     
if [ -f $FILE ]; then
   echo "File $FILE exists, cleaning"
   make Clean
else
   echo "File $FILE does not exist, no cleaning needed"
fi

if [ "$MACHINE" == "Linux" ];
then
	#BUILD_FILE='linux_amd64_pgf77_pace'
	module load mvapich2/2.3.7-1 
	#module load netcdf-fortran/4.6.1-mva2-hdf5-1.14
	export MPI_HOME='/usr/local/pace-apps/spack/packages/linux-rhel9-x86_64_v3/gcc-12.3.0/mvapich2-2.3.7-1-qv3gjagtbx5e3rlbdy6iy2sfczryftyt/'
	#export NETCDF_HOME='/usr/local/pace-apps/spack/packages/linux-rhel9-x86_64_v3/gcc-12.3.0/netcdf-c-4.9.2-hv6rtvb7476ibpgxpi54pcgkbouzswfl/'
	BUILD_FILE='linux_amd64_gfortran'
else
	BUILD_FILE='darwin_amd64_gfortran'
fi

if [ -z $OPT2 ]; then
   echo "===== Building with ${BUILD_FILE} at $(pwd) ====="
   $ROOT/tools/genmake2 -mods ../code -optfile $ROOT/tools/build_options/$BUILD_FILE -rootdir $ROOT
elif [ "$2" == "-mpi" ]; then
   echo "====== Building with MPI and ${BUILD_FILE} at $(pwd) ====="
   $ROOT/tools/genmake2 -mods ../code -mpi -optfile $ROOT/tools/build_options/$BUILD_FILE -rootdir $ROOT
else
   echo "$2"
   echo "invalid 2nd argument"
   exit 2
fi

echo " Done compiling, moving to make depend..."
make depend -s
echo " Done with make depend, moving to make..."
make -sj 4

if [ "$MACHINE" == "Mac" ];
then
   afplay /System/Library/Sounds/Funk.aiff &
fi 

echo "(⌐■_■) Done building"
