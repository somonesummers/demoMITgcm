To install packages from Paul's scripts run:

You must install conda (https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

Then open a terminal, navigate to this directory (where this readme is) then run:


conda env create -f environment.yml

press 'y' and you should be good to go!


(if you want to save your new env, you can use 
	conda env export | grep -v "^prefix: " > environment.yml
)