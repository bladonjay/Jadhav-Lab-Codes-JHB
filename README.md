# Jadhav Lab Codes
 Code repositories from Jadhav Lab
 This repository has all my code from the jadhav lab.  Right now there are three projects
 
  1. BetaOdorProject: olfactory guided T maze task, this was published under DOI: https://doi.org/10.7554/eLife.79545
     This repository also contains an archived version of 'general_repo' which will need to be in your matlab search path.
     To start analysis run the function: Jadhav-Lab-Codes\BetaOdorProject\figures1-6Functions\PipelineScriptsSTARTHERE\cs_Pipeline
     To start analysis for figure 6: cd C:\Users\Jadhavlab\Documents\gitRepos\Jadhav-Lab-Codes\BetaOdorProject\figure7Functions
     then run 'A_Figure7Wrapper'
     
     The datafiles to run this code live here: 
     https://figshare.com/articles/dataset/Rhythmic_coordination_of_hippocampal-prefrontal_ensembles_for_odor-place_associative_memory_and_decision_making/19620783
     ** the LFP data for the last animal is too large for the 20 gb limit in Figshare.  You can email me at jhbladon(at)brandeis(dot)edu for those data or wait
     until they are uploaded to the dandi archive.
    
     For the 'pipelineScripts' functions you will need to change some directory-finding functions, namely JHB_setBetaOdorpaths and cs_setPaths
     For the 'Wrapper' function, just bring the 'BetaOdorData' file (~6gb) into your workspace and make sure you have in your path the folder called general-repo
     
     
  
  2. FMR1 project: This code is not fully developed, but it has scripting for deeplabcut, SpikeGadgets data preprocessing and
      analysis of two behavioral paradigms: the social W and the W-maze home switch paradigm. Not for public use yet
      
  3. DG-GABAa5 KO animals: This codebase is a partial workup of some data collected by Elif Engin at Harvard. This is an abbreviated version of the radial
     arm maze, wherein there are 4 outbound arms and a single home arm.  The single out arm maze is shifted each time.
      
 
 Be forewarned these code bases often have dependencies and query hardcoded filepaths.  The second two items (FMR1 and DG-GABAa5) are not for public use.
