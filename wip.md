# Virome pipeline WIP
This repository is a work-in-progress where everyone can share their scripts of the pipeline in their own folder and where we can work together on one 'consensus' pipeline.

To do this everyone can make their own branch of this repo, change in your branched consensus script what you think could be better and then create a pull request to merge these changes in the main branch. Subsequently, everyone can review the changes that will happen on the merge and agree/disagree. When everyone agrees the pull request can be merged in the main branch. This way we have a clean way to build a 'master' pipeline from all versions that are currently circulating in the lab.

A ppt with introduction to Github is also included in this repo, I recommend everyone to go over it before working on this project. Be aware that Github changed its 'vocabulary', the 'main' branch used to be 'master' branch. This isn't updated in the ppt.

This repo also contains a yaml file from which you can build your conda environment that will be used to run the pipeline. If you haven't installed Miniconda yet, you can find a guide [here](https://github.com/Matthijnssenslab/virome-pipeline-wip/blob/main/installing-miniconda/installing-miniconda.md).

## Setup
```bash
git clone https://github.com/Matthijnssenslab/virome-pipeline-wip.git
cd virome-pipeline-wip/old-versions
mkdir <your-name> #copy your old version(s) of the pipeline here
git add .
git commit -m "versions <your-name>"
git push
git checkout -b <your-name> #will create a new branch and you will move into this new branch
git push --set-upstream origin <your-name> #will push your branch to the remote repository
```
## Setup conda environment
```bash
cd virome-pipeline-wip
git checkout Lander
conda env create -f viper.yml
#Install viper.sh and Cluster_genomes.pl in conda env PATH
cd bin
ln -sfr * $CONDA_PREFIX/envs/viper/bin
#link Krona databases
rm -rf $CONDA_PREFIX/envs/viper/opt/krona/taxonomy
ln -sf /staging/leuven/stg_00029/DB/Krona/latest/ $CONDA_PREFIX/envs/viper/opt/krona/taxonomy
conda activate viper
```
## Update local repo
Keep your local repository up to date by **often** pulling the repository on Github!
```bash
cd virome-pipeline-wip
git pull
```
## Making and committing changes
Always make sure that you are working in your own branch!
```bash
git checkout <your-name> #will move to your branch
git status #will tell on which branch you are
```
Just like you have to update your local repo with `git pull` often, also commit your changes **often** to avoid compatibility errors with the remote repository.
```bash
git add . #this will add all changed files to the staging area of git
#or
git add <specific-file> #this will only add the specific file
git commit -m "meaningful message" #will commit your changes, add a meaningful message to keep better track of your changes
git push #will push the changes to the remote repository
```
## Overview of git commands
<img src="https://material.bits.vib.be/topics/git-introduction/images/conceptual_areas_push.png"/>
