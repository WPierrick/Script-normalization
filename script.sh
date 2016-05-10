#!/bin/sh

qsub -b y -cwd -l vf=150G -l h_vmem=150G Rscript RUV_Script_03.R

#  Script.sh
#  
#
#  Created by Pierrick Wainschtein on 27/04/2016.
#
