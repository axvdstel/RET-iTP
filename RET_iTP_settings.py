#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 09:58:53 2021

@author: anne-xandervanderstel
"""

#########    USER PARAMETERS     ##################

# 3' Adapter fixed sequence to be used for adapter trimming
adapter_sequence = ''
# number of ambigeous bases of 3' adapter
adapter_extra = 6

#number of nucleotides protected by the ribosome
ribosome_protect = 17

# gives file names of different RAW NGS files placed in the 'RAW_reads' folder, please unzip these files
# Note! give the path to the reverse reads file

RET_1_path = 'RAW_reads/'
RET_2_path = 'RAW_reads/'
control_1_path = 'RAW_reads/'
control_2_path = 'RAW_reads/'
library_path = 'RAW_reads/'


#download and put into the correct directories the fasta genome and GFF file of your organism
genome_file = 'genome_files/'
GFF_file = 'genome_files/'


# The following executables need to be installed and their respectives paths need to be gives here

# minimap2 (https://github.com/lh3/minimap2)
minimap2_path = ''
#BBDUK from BBmap(https://sourceforge.net/projects/bbmap/)
BBDUK_path = ''



#computing_cores_tobe_used
cores = 4

