#!/usr/bin/env bash
# Script for running MITgcmuv
# Paul Summers Sept 2025
# invoked as:
# bash makeRun.sh [-mpi]
# -mpi for mpi running. This requires specific packages to be loaded, currently this is set up for GaTech Pace
# MPI_HOME will need to be adjusted for local configuration. 
# if you list some
set -e
OPT1=$1
OPT2=$2
TIME="$(date +"%y%m%d_%H%M%S")" 
TIMENICE="$(date +"%H:%M:%S")" 
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     MACHINE="Linux";;
    Darwin*)    MACHINE="Mac";;
esac
echo "Idenitfied machine as ${MACHINE}"
echo "Already built: clean up run folder, then make simlinks and run"

if [ -d results ]; then
	cd results
else
	echo "results folder missing, making now"
	mkdir results
	cd results
fi

touch test.txt #this ensures the dir is not empty (rm -f * also probably works, but hey if it ain't broke, don't fix it)
rm *
ln -s ../input/* .
cp ../build/mitgcmuv .

echo "Ready to run -=三(ง ˙o˙)ว"

echo "Running from  $(pwd) at $TIMENICE"
if [ "$MACHINE" == "Mac" ];
then
	if [ -z $OPT1 ]; then #standard run
	   time ./mitgcmuv > ../Report$TIME.txt
	elif [ $OPT1 == "-mpi" ]; then
		if [ -z $OPT2 ]; then
	  		echo "not enough arguments, need to specify number of cores"
	  		exit 2
  		else
  			time mpirun -np $OPT2 ./mitgcmuv > Report$TIME.txt
  		fi
	fi
	cd ..
	if [ -d figs ]; then
   		echo "figs folder already made"
   	else
   		echo "figs folder missing, making now"
   		mkdir figs
   	fi	
	afplay /System/Library/Sounds/Funk.aiff
else
	if [ -z $OPT1 ]; then #standard run
		./mitgcmuv
	elif [ $OPT1 == "-mpi" ]; then
		srun ./mitgcmuv
	else
		echo "invalid argument: $OPT1"
		exit 2
	fi
fi

echo "Done running  ദി(˵•̀ᴗ-˵)✧"
