Paul Summers, September 2025. A demo of MITgcm with iceberg2 and iceplume for group meeting.
Dependencies loaded with conda. 
Validated on Mac Silicon (M series chips). 
Email paul.summers@gatech.edu with any questions

To install packages from Paul's scripts run:

You must install conda or miniconda (https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

	Then open a terminal, navigate to this directory (where this readme is) then run:
	
conda env create -f environment.yml

	This will install a new conda env, if it fails try to resolve these things, but hopefully it'll say it was successful. It this fails, please see my tips at the very bottom

(if you want to re-save your new env, you can use 
	conda env export | grep -v "^prefix: " > environment.yml
)

conda activate demo
	
	This will activate the new conda env which is called 'demo'. This has everything you need to run my python scripts, as well as even the correct version of gfortran. I'm not convinced this will work for everyone, but so far it has worked on all Mac Silicon machines I have tried. If you're doing this on a cluster, they should have gfortran already, talk with your support team if not.

cd GenerationScripts

	This is moves you to where the generator script lives

python finModel.py

	This is how you run the generator script. It will ask you to confirm the description. It will then produce numerous plots confirming to you what is being created. You must close these plots to keep going.  

cd ../experiments/bergDemo

	This moves to the newly created folder. All of the following commands are run from the experiment folder. You will never have to move terminal window to another folder for building/running/visualizing. 

bash makeBuild.sh
bash makeRun.sh

	This builds, (which will return many warnings and comments) then runs (which also has a few warnings). When done with building and/or running you'll hear some funk sounds and get a happy ascii art face (aka (⌐■_■) ). 
	Running should report "STOP NORMAL END" when done, as well as some stats on how long it took to run. 

python ../giflookSideAvg.py
python ../giflookMap.py -z -100 -s
python ../giflookXSlice.py -x 9000

	This will run plotting, saving figures to the figs/ folder. The first plots the side view of the fjord, averaged across the fjord. The second plots a map view at depth of -100 meters, with a shadow for ice bergs (-s for shadow). The third plots an across fjord slide at 9000 meters (-x 9000) away from the glacier. You can add "-q" to only save pngs of the final timestep or "-qq" to only save pngs and view the images in the matplotlib native figure viewer. "-h" will show you the full options for these files. They load some defaults from files in "input/plotHelper.py" like the colorbar values. You can change those within that file (input/plotHelper.py) if you'd like.

If you find this fun, you can change the input/data file to have an endTime of 864000 (10 days in seconds). Next you will call makeRun.sh againe (no need to makeBuild again) and it should take around 5 minutes to run. In this 10 day run you can see the 2 separate outflows from the plume/bergs separating. You will have to make the figures again to see this. 'python ../giflookXSlice.py -x 9000' shows this most clearly. 
	

	
IF CONDA FAILED WHEN LOADING MY ENV:
	Here are instructions for manually making a similar env. 
	You should make a new blank conda env, and then activate it

conda create --name demo
	
	You may have to reply 'y' here to continue

conda activate demo

	Then installing these packages via command line:

conda install conda-forge::gfortran anaconda::python conda-forge::matplotlib conda-forge::imagemagick conda-forge::gsw conda-forge::cmocean conda-forge::scipy conda-forge::netcdf4

	again you may have to reply 'y' to continue	

pip install MITgcmutils

	Once those are all successful I believe you can now run everything! 
