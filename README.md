# VASP FILE ANALYSIS

Code written by C. Barger - 2018

input : VASP FILES = OSZICAR, CONTCAR and EIGENVAL
	Additional file needed will be NAME.txt, which should be created in the shell script.
	
	A directory called "Crystal" exists within the directory for the code. Inside of this directory should be the files necessary to analyze the crystal model's data using the code. This include the crystal's CONTCAR, EIGENVAL, and OSZICAR output files from VASP. In addition, the code generates the IS_CRYSTAL.txt file to put inside of the "crystal" directory if it does not exist already. This file is necessary so that the code doesn't treat the crystal model as having hydrogen atoms.
	
	In the RUN executable, there is an example that shows you how to move around the files into the Nano-Pore-Analyzer and run the code for a subset of models. This subset is commented out. Also in the RUN file there is a section that executes the Nano-Pore-Analyzer code for the crystal model. There should be little change, if any change at all, to this section of the RUN executable.

to use: Compile the code using the make command
   	The code is run using ./hello
	The code should write information to DATA.csv

output :  DATA.csv which includes columns for features = Band Gap, atomic density ratio values, atomic densities, etc

script: RUN file allows one to run code for multiple systems and place results in one DATA.csv file. The RUN file should be modified by 	the user such that the code is executed for all of the needed models.

