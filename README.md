# Jadhav Lab Codes
 Code repositories from Jadhav Lab
 This repository has all my code from the jadhav lab.  Right now there are three projects
 
  1. BetaOdorProject: olfactory guided T maze task, this was published under DOI: https://doi.org/10.7554/eLife.79545
     This repository will require also the download of 'general_repo' and it needs to be in your matlab search path.
     To start analysis run the function: Jadhav-Lab-Codes\BetaOdorProject\figures1-6Functions\PipelineScriptsSTARTHERE
     To start analysis for figure 6 run: cd C:\Users\Jadhavlab\Documents\gitRepos\Jadhav-Lab-Codes\BetaOdorProject\figure7Functions
     then run 'A_Figure7Wrapper'
     
     The datafiles to run this code live here: 
     https://figshare.com/articles/dataset/Rhythmic_coordination_of_hippocampal-prefrontal_ensembles_for_odor-place_associative_memory_and_decision_making/19620783
     For the 'pipelineScripts' function you will need to change some directory-finding functions, put error checking on
     For the 'Wrapper' function, just bring the 'BetaOdorData' file (~6gb) into your workspace
     
     
  
  2. FMR1 project: This code is not fully developed, but it has scripting for deeplabcut, SpikeGadgets data preprocessing and
      analysis of two behavioral paradigms: the social W and the W-maze home switch paradigm. Not for public use yet
      
  3. DG-GABAa5 KO animals: this is an abbreviated version of the radial arm maze, wherein there are 4 outbound arms and a single
      home arm.  The single out arm maze is shifted each time.
      
 
 These code bases often have dependencies and query hardcoded filepaths.  They are not intended for general use by any means.
