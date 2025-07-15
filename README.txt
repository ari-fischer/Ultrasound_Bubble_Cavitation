15 July 2025

Code and Data Repository Accompanying Ultrasound-Induced Bubble Cavitation Simulations

This directory contains custom codes written to simulate the cavitation of Ar microbubbles in water exposed to ultrasound irradiation. These codes were written by Ari F. Fischer, Robert C. Maligon Querimit, and Tej S. Choksi, and implement the model originally published in 10.1098/rspa.2001.0784. These codes were used to generate the results reported in a manuscript submitted from the research group of Asst. Prof. Tej S. Choksi at Nanyang Technological University in Singapore.

All codes are implemented in MATLAB, and were either run (i) locally on a Macbook air or (ii) remotely on using high-performance computing (HPC) resources. 

--The "custom_functions" directory contains custom MATLAB functions that run the bubble simulations and process the data as used in the submitted manuscript. 

--The "NASA_Polynomials" directory contains the coefficients for NASA polynomials used to calculate thermodynamic properties of the gas components present in bubbles during collapse

--The "Reaction_Parameters" directory provides details of the reaction network considered for H2O consumption. These include the rate constants for each elementary reaction, the extents of reactions for those steps, and the stoichiometric coefficients for the forward and reverse rates of each step.

--The "run_analysis_scripts" executes the codes in the "custom_functions" directory in order to generate the datasets reported in the submitted manuscript. The large data-sets requiring HPC are reported in the "figures_main_HPC" directory, and smaller data-sets collected locally are reported in the "figures_main_local" directory. The codes and data are organised by the figure number in which they are reported in the submitted manuscript. 

For codes run locally, simulation results are reported as MATLAB variables in ".mat" files, readable within MATLAB--some data are also reported in ".csv" files. 

For codes run with HPC, datasets are reported as ".csv" files. 

The MATLAB codes used to generate these figures are also included. For figures plotting data from ".mat" files, these can be reproduced by first loading the ".mat" variables into the MATLAB workspace, then running the figure-plotting script. The figures plotting data from ".csv" files can read the data provided a correct path to the data.


