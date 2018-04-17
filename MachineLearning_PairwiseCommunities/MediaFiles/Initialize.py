# 17 January 2017
# Author: Demetrius DiMucci

#!/usr/bin/python

# This script does the following:

# 1. Identify all of the models in the target directory and put them into a list.
#       find_mod() does this.
# 2. Create an unique output folder for each experiment
#       naming() does this.
# 3. Populate the folder with the required files.
#       move_files() does this.

# Native modules
import sys #for reading from the command line
import os #for making regular bash calls
import time #for timing the runs to get a feeling for how long everything will take
import subprocess
from subprocess import call #for making a bunch of threads to run new COMETS instances simultaneously
import itertools
import glob

## 
## MODIFY THESE PATHS TO POINT TO WHERE YOUR MODELS LIVE
## I did this on my institution's cluster
Models_path = '/projectnb/cometsfba/dimucci/Pairwise/Models'
script_path = '/projectnb/cometsfba/dimucci/Pairwise/Scripts'
media_path = '/projectnb/cometsfba/dimucci/Pairwise/MediaFiles'
path = '/projectnb/cometsfba/dimucci/Pairwise/Output_301/'
#path = './'
##

# Tracking variable
folderCount = 0
moveCount = 0
# Define the functions

def find_mod(path):
        #this function finds all the model files in the given path (probably current) and returns them so they can all be run through COMETS
        flist = glob.glob( path +'/*xml.mat') #glob-glob anything with the quoted extension. These are our model files. 
        mod_list = []   # a list to store them in. we need to do some parsing before they're in the form we need.
        #print len(flist)
        for mod in flist: #for each model
                mod = mod.strip().split('/') # get the part after the "/"
                mod = mod[len(mod)-1]
                mod_list.append(mod) # put that part in our real models to run list.
        return mod_list #return this list.

# This function will create a list of folder names make the directory and copy the models into it.
def nameFolders(allCombos,Models_path):
        folder_list = []
        for combo in allCombos:
                models  =  sorted(combo)
                folder_name = path
                for model in models:
                        model = model.split('_xml.mat')[0]
                        folder_name += '_' + model
                folder_list.append(folder_name)
                os.system('mkdir ' + folder_name)
                for model in models:
                        os.system('cp ' + Models_path + '/' + model + ' ' + folder_name)
        return folder_list

def moveFiles(folder_list, path):
        moveCount =0
        for folder in folder_list:
                os.system('cp '+ path + '/comets_scr'+' '+folder)
                os.system('cp '+ path + '/comets_script'+' '+folder)
                os.system('cp '+ path + '/makeLayouts.m'+' '+folder)
                os.system('cp '+ path + '/global_params.txt'+' '+folder)
                os.system('cp '+ path + '/package_params.txt'+' '+folder)
                os.system('cp '+ path + '/runScript'+' '+folder)
                os.system('cp '+ path + '/makeLayouts.r'+' '+folder)
                os.system('cp '+ path + '/compounds.m'+' '+folder)
                os.system('cp '+ path + '/ratios.m'+' '+folder)
                moveCount += 1
        return moveCount


  #          #          #
   #         #         #
    #        #        #
     #################
#####      MAIN       #####
     #################
    #        #        #
   #         #         #
  #          #          #


mod_list = find_mod(Models_path)

# # This block will generate a list of all model combinations
bothCombos = []
allCombos = []
combos = tuple(itertools.combinations(mod_list,2))
monoCombos = tuple(itertools.combinations(mod_list,1))
bothCombos.append(monoCombos)
bothCombos.append(combos)

for i in range(0,len(bothCombos)):
        for j in range(0,len(bothCombos[i])):
                allCombos.append(bothCombos[i][j])

folder_list = nameFolders(allCombos,Models_path)
moveFiles(folder_list,media_path)

