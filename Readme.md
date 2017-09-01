# CHMD
Automatic reconstruction of the muscle fibres from the superficial layer fibres data

## Getting the source code
CHMD has two remote repositories: one private at https://forge.kiv.zcu.cz:3443/git/cholt-cadaver-fibres.git and one public at https://github.com/kivzcu/CHMD.git. The private repository contains two main branches, master and CHMD_master, whereas CHMD_master contains the public stuff and master contains the proprietary code. The public repository contains one main branch named master that corresponds to the branch CHMD_master of the private repository. If you are neither regular or associate member of the research team at the School of Computing, University of West Bohemia, you have only access to the public repository, otherwise, it is recommended to work with both remote repositories but one local repository. How to do it?

1) Using Git Extensions clone the private repository (master branch)
2) Create a new branch at the latest commit of origin/CHMD_master with name CHMD_master and check it out
3) Add a new remote repository named GitHub and Url https://github.com/kivzcu/CHMD.git
4) Pull the changes from the GitHub (and push them to the private repository origin/CHMD_master)
5) Checkout master branch (the private one)

6) Cherry-pick commits from the master that should be made public for CHMD_master and cherry-pick public commits that are important for the private master branch. 

