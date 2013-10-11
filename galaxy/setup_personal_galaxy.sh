#!/bin/bash
#
# Check out a new personal galaxy 
#
# Based on https://projects.uabgrid.uab.edu/galaxy/wiki/GalaxyDevelopment
#

PROJECT_DIR=$HOME/projects/galaxy
GALAXY_DIR=$PROJECT_DIR/galaxy
if [ -e  $GALAXY_DIR ]; then
    echo "$GALAXY_DIR already exists; remove it before running a clean checkout"
    exit 1
fi

# load needed modules 
source /etc/profile.d/modules.sh
echo module load galaxy/galaxy-developer
module load galaxy/galaxy-developer

mkdir -p $PROJECT_DIR
cd $PROJECT_DIR

echo git clone ssh://${USER}@git.uabgrid.uab.edu/home/git/repositories/galaxy.git
git clone ssh://${USER}@git.uabgrid.uab.edu/home/git/repositories/galaxy.git

git remote add rollout ssh://galaxy@galaxy.uabgrid.uab.edu/share/apps/galaxy/RollOutRepo/galaxy.git

echo cd $HOME/projects/galaxy/galaxy
cd $HOME/projects/galaxy/galaxy
echo git checkout develop
git checkout develop

echo "Suggested port: "`grep $USER /etc/passwd | cut -d : -f 3`" (user id for $USER)"
echo ./personalize.sh 
./personalize.sh 

echo ./run.sh --daemon
./run.sh --daemon

