Mountainlab-JS for JadhavLab
======

Repository for Jadhav Lab code for the use of MountainSort-JS

TODO: convert franklab script to convert mountainsort output to FF files (spikes & cellinfo)
TODO: launcher script to simplify qt-mountainview launch

Automatic Setup (9-24-18)
------
* Clone this repo
* In a terminal, cd into this repo
* run  setup_mountainlab.sh
    * `. setup_mountainlab.sh`
* open  `~/conda/env/mlab/etc/mountainlab/mountainlab.env` in text editor of your choice
* add `ML_TEMPORARY_DIRECTORY='/path/to/data/drive/tmp/mountainlab-tmp'`
* To your matlab path, add `~/conda/env/mlab/lib/node_modules/mountainlab/utilities/matlab/mdaio/`


Manual Setup
------
* Download and install miniconda

```shell
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
        bash miniconda3.sh -bp ~/conda
        echo ". ~/conda/etc/profile.d/conda.sh" >> ~/.bashrc
```

* Setup a conda environment and install mountainlab and all processors

```shell
        conda create --name mlab
        conda activate mlab
        conda update conda
        conda config --set max_shlvl 1
        conda install -c flatiron -c conda-forge mountainlab mountainlab_pytools ml_ephys ml_ms4alg ml_ms3 ml_pyms qt-mountainview
```
* Add `~/conda/env/mlab/lib/node_modules/mountainlab/utilities/matlab/mdaio/` to your matlab path

* Get additional processors from franklab for drift tracking and tagged curation
    * Actually copy or symlink the folders from this repository into `~/.mountainlab/packages/`
    * [FrankLab Tagged Curation](https://bitbucket.org/franklab/franklab_mstaggedcuration/src/master/)
    * [FrankLab Drift Tracking](https://bitbucket.org/franklab/franklab_msdrift/src/master/)
    * These should be cloned into `~/.mountainlab/packages/`

* Create configuration file for mountainlab
    * `touch ~/conda/env/mlab/etc/mountainlab/mountainlab.env`
    * Now you can set the temporary directory path to be on the same drive as your data
    * Also add the `~/.mountainlab/packages/` to the mountainlab package search path
    * Just modify and add these lines to the `mountainlab.env` file you created:
        * `ML_TEMPORARY_DIRECTORY='/path/to/data/drive/tmp/mountainlab-tmp'`
        * `ML_ADDITIONAL_PACKAGE_SEARCH_DIRECTORIES='~/.mountainlab/packages'`

### Updates 9-23-18
The franklab msdrift and mstaggedcuration packages will throw errors as is due to changes in package locations in the new mountainlab-js. So instead copy the `franklab_mstaggedcuration` and `franklab_msdrift` folder from this repository into `~/.mountainlab/packages/`. 

Usage
------
* run exportmda on your raw trodes data
* Required Directory Structure (+ - required naming scheme)
```
    animal_container_folder
    --> raw_data_folder (eg. RZ9)
        --> day_directories +(day#_date eg. 01_180924)
            --> Raw data files, etracted binaries, etc. +(animID_day#_date_whatever eg. RZ9_01_180924_1Sleep)
    --> direct_folder +(animID_direct eg. RZ9_direct) 
```
    * required naming scheme is because currently the scripts parse animID, day #, date and tetrode number from file names
* ml_process_animal(animID,rawDir) will mountainsort all days for that animal. Output will be saved to animID_direct/MountainSort in folders for each day, with subfolders for each tetrode
    * I recommend passing this function a tet_list to restrict the tetrodes you run it on to only those you want to cluster, since MountainSort will take ~3.5GB per tetrode for 1hr of recording
* This wrapper function will:
    * create the output directories
    * make prv links to the raw mda files and copy over the timestamps.mda (TODO: make it use a link to save space)
    * Bandpass filter (600-6000 Hz), mask artifacts and whiten  the raw data. Saved as filt.mda (filtered and masked), and pre.mda (whitened filt.mda).
    * Split into epochs and spike sort each epoch separately
    * Track drift and combine clusters across epochs (output: firings_raw.mda)
    * Calculate cluster and isolation metrics (output: metrics_raw.json)
    * Automatically tag clusters to aid with manual verification (output: metrics_tagged.json), threshold taken from frank lab
* You can manually run each step and customize parameters (i.e. filtering band, thresholds, etc) using ml_filt_mask_whiten and ml_sort_on_segs which take the tetrode results directory as the input (animID_direct/MountainSort/day_folder.mountain/tetrode_folder)
    * read the help for each function to figure out how to customize variables in the function call
* additionally, all parameters  can be overriden using the params.json file inside the tetrode results directory by including fields in json format.

qt-mountainview
------
* cd into the tetrode directory you want to view
* launch qt-mountainview
    `qt-mountainview --raw=raw.mda.prv --filt=filt.mda --pre=pre.mda --firings=firings_raw.mda --cluster_metrics=metrics_tagged.json --samplerate=30000`
