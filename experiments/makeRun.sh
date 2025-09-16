#!/usr/bin/env bash
set -e
TIME="$(date +"%y%m%d_%H%M%S")" 
TIMENICE="$(date +"%H:%M:%S")" 
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     MACHINE="Linux";;
    Darwin*)    MACHINE="Mac";;
esac
echo "Idenitfied machine as ${MACHINE}"
echo "Already built: clean up run folder, then make simlinks and run"
cd results
touch test.txt #this ensures the dir is not empty (rm -f * also probably works, but hey if it ain't broke, don't fix it)
rm *
ln -s ../input/* .
cp ../build/mitgcmuv .

echo "Ready to run -=三(ง ˙o˙)ว"

echo "Running from  $(pwd) at $TIMENICE"
if [ "$MACHINE" == "Mac" ];
then
	time ./mitgcmuv > ../Report$TIME.txt
	cd ..
	if [ -d figs ]; then
   		echo "figs folder already made"
   	else
   		echo "figs folder missing, making now"
   		mkdir figs
   	fi	
   	for NAME in "$@"
		do
			python $NAME
		done	
	afplay /System/Library/Sounds/Funk.aiff &
else
	./mitgcmuv
fi

echo "Done running  ദി(˵•̀ᴗ-˵)✧"
