#!/bin/bash

if ! [ -x "S(command -v conda)" ]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/Downloads/miniconda3.sh && \
    sudo chmod 777 ~/Downloads/miniconda3.sh && \
    bash ~/Downloads/miniconda3.sh -bp ~/conda && \
    echo ". $HOME/conda/etc/profile.d/conda.sh" >> ~/.bashrc
    . ~/.bashrc
fi
echo "Setting up mlab conda environment..."
conda create --name mlab && \
conda activate mlab && \
conda config --set max_shlvl 1 && \
echo "Installing mountainlab, qt-mountainview and associated processors..."  && \
conda install -c flatiron -c conda-forge mountainlab mountainlab_pytools ml_ephys ml_pyms ml_ms3 ml_ms4alg qt-mountainview && \
echo "Setting up mountainlab environment..." && \
local mountainENV=~/conda/envs/mlab/etc/mountainlab/mountainlab.env 
touch $mountainENV
echo "ML_ADDITIONAL_PACKAGE_SEARCH_DIRECTORIES='~/.mountainlab/packages'" >> $mountainENV &&  \
ml-config
echo "Copying franklab_msdrift and franklab_mstaggedcuration to ~/.mountainlab/packages/"
mkdir ~/.mountainlab/packages
cp -R franklab_msdrift ~/.mountainlab/packages/franklab_msdrift
cp -R franklab_mstaggedcuration ~/.mountainlab/packages/franklab_mstaggedcuration
echo "set ML_TEMPORARY_DIRECTORY in $mountainENV to be on the same drive as the data for ease of processing and space management" 
if [ -x "$(command -v lolcat)" ] && [ -x "$(command -v figlet)" ]; then
    echo 'mountainlab setup complete' | figlet | lolcat
    echo 'add ~/conda/envs/mlab/lib/node_modules/mountainlab/utilities/matlab/mdaio/ to matlab path' | lolcat
else
    echo 'mountainlab setup complete'
    echo 'add ~/conda/envs/mlab/lib/node_modules/mountainlab/utilities/matlab/mdaio/ to matlab path'
fi

