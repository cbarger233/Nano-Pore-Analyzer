# VASP FILE ANALYSIS

Code written by C. Barger - 2018

input : VASP FILES = OSZICAR, CONTCAR and EIGENVAL
	Additional file needed will be NAME.txt, which should be created in the shell script.

to use: Compile the code using the make command
   	The code is run using ./pore
	The code should write information to DATA.csv

output :  DATA.csv which includes columns for features = Band Gap, atomic density ratio values, atomic densities, etc

script: RUN file allows one to run code for multiple systems and place results in one DATA.csv file

