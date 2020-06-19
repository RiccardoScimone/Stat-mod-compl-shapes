This code allows to perform a complete simulation step as the ones presented in the paper. 
Point cloud handling is managed resorting to CloudCompare (http://www.danielgm.net/cc/release/)
The following preliminary steps are required to run correctly the code.

The directory CloudCompareBinaries contains the CloudCompare executable for Windows 64 bits, so it should already work correctly on 
this OS. If you have a different OS, these are the steps to follow

(1) Download the compressed archive 7z version of CloudCompare, compatible with your OS, from http://www.danielgm.net/cc/release/
(2) Clear the content of CloudCompareBinaries directory and unzip the downloaded archive there. At the end, CloudCompareBinaries should con
tain a CloudCompare executable.

Now the directory contains a CloudCompare executable which can be called from R through command-line mode.

The code is organized as follows:

(1) The directory MatlabScriptEgg2D contains matlab code, which is meant to be called from R,
 to generate the shape point cloud. The user is supposed to modify this scripts in order to set the simulation parameters.
(2) The directory RscriptSfpca contains R code to call matlab and generate shapes, and to perform all the analysis proposed
 in the work.

We give an example on how to proceed for a simulation run. 
Suppose we choose to simulate for Scenario I (II, III). This are the step to be followed.

(1) Go to the MatlabScriptEgg2D directory. Open the run_from_r.m script e uncomment the line
 Scenario_I_Simulation (Scenario_II_Simulation, Scenario_III_Simulation), commenting the 
others if they are not commented.

(2) Open the simulation script Scenario_I_Simulation.m (Scenario_II_Simulation.m, Scenario_III_Simulation.m)
 and change the corresponding parameter mu_s1 (eta, nu) as desired. This parameters are always set at the
beggining of the scripts, followed by the comment %TO BE SET.

(3) Go to the RscriptSfpca directory, open distance_computation.r, and change the parameter
choose_Scenario with the number of the chosen scenario (I, II, III)
(4) Open Chart_generator_pca.r and do exactly the same thing
(5) If not already installed, install the r packages foreach, doParallel, fda, svMisc, matlabr
(6) You are now ready, open run.r and execute all the lines. The images resulting from the analysis
will be saved in the corrisponding img directory.



If you want to visualize the Point Clouds, you can uncomment the corresponding lines, making plots, in the matlab scripts. 