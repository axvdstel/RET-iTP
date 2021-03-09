#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 14:25:27 2021

@author: anne-xandervanderstel
"""
import os
import subprocess
import Bio
import multiprocessing
import defaultdict
import re
from Bio import SeqIO
import numpy as np
from RET_iTP_functions import *
from RET_iTP_settings import *
import pandas as pd
#########  RET-iTP Main script  #################


def log(string):
    with open('logfile.txt', 'a') as f:
        f.write(str(string) + '\n')
    print(str(string))
    return




#######################  INITITALIZING   #################################################
log(Start Initializitaion...)

sample_list = [RET_1_path, RET_2_path, control_1_path, control_2_path, library_path]

#looking for innput files
for file in sample_list:
    if os.path.exists("%s" %(file)) == False:
        print('uhuhm, no raw files detected!! Put the correct files in the correct folder and adjust the script!')
    if file.split('.')[-1] != 'fq':
        print('Please unzip the RAW files')

for file in [RET_1_path, GFF_file]:
    if os.path.exists("%s" %(file)) == False:
        print('uhuhm, no genome files detected!! Put the correct files in the correct folder and adjust the script!')

os.mkdir('Trimmed')
os.mkdir('Mapped')
os.mkdir('Analysis')


def genome():

    for genome in SeqIO.parse(genome_file, 'fasta'):
        return genome

def CDS():
    with open(GFF_file, 'r') as f:
        features_plus = {}
        features_min = {}
        for x in f.read().split('\n')[:-1]:
            if x[0] != '#':
                if x.split('\t')[2] == 'CDS':
                    if x.split('\t')[6] == '+':
                        features_plus[(int(x.split('\t')[3]),int(x.split('\t')[4]))] = ([y for y in x.split('\t')])
                    elif x.split('\t')[6] == '-':
                        features_min[(int(x.split('\t')[3]),int(x.split('\t')[4]))] = ([y for y in x.split('\t')])
    return features_plus, features_min

genome = genome()
genome_l=len(genome.seq)
CDS_plus, CDS_min = CDS()




log(Initiatlization complete)

#################### ADAPTER CHOPPPING ##################################
log('STARTING ADAPTER_TRIM')


for sample_path in sample_list:
    rev_trimmer(sample_path)

log('Trimming finished')

########### MAPPING reads to genome




for sample_path in sample_list:
    minimapper(sample_path)

log('Mapping finished')


########## Make WIG files ###############




pool = multiprocessing.Pool(processes=cores)
out = pool.map(find_Psites, sample_list)
pool.close()
WIG_writer_background()
log('WIGs writing finished')

####### Smooth library reads

bg_plus, bg_min =background_parser()
bg_plus_smooth=smoother(bg_plus)
bg_min_smooth=smoother(bg_min)
WIG_writer('Analysis/background_plus_smooth.wig', bg_plus_smooth, genome)
WIG_writer('Analysis/background_min_smooth.wig', bg_min_smooth, genome)

log('Library smoothing finished')


########### Correct reads for library



bg_freq_plus = frequency_reads(bg_plus_smooth)
bg_freq_min = frequency_reads(bg_min_smooth)
bg_freq_plus = bg_lim(data_new)
bg_freq_min = bg_lim(data_new)

for sample_path in sample_list:
    if sample_path == library_path:
        continue
    sample_name = sample_path.split('/')[-1]
    handle_plus = 'Analysis/RET-iTP_reads_%s_plus.wig' %(sample_name)
    handle_min = 'Analysis/RET-iTP_reads_%s_min.wig' %(sample_name)
    reads_plus = read_WIG(handle_plus)
    reads_min = read_WIG(handle_min)

    scores_plus = frequency_reads(reads_plus) / bg_freq_plus
    scores_min = frequency_reads(reads_min) / bg_freq_min
    WIG_writer('Analysis/RET-iTP_scores_%s_plus.wig' %(sample_name), scores_plus, genome)
    WIG_writer('Analysis/RET_iTP_scores_%s_min.wig' %(sample_name), scores_min, genome)

log('Library coorection finished')
################# Score the data on nucleotide level


all_scores = RET_iTP_scorer('Analysis/RET-iTP_scores_all.csv')

log('scoring finished')
############### FIND ALL POSSIBLE ORFs


pool = multiprocessing.Pool(processes=2)
orfs_plus_coordinates, orfs_min_coordinates_rc = pool.map(orf_finder, [genome.seq, revcom(genome.seq)])
pool.close()
orfs_min_coordinates = defaultdict(lambda:1,{genome_l-x:[genome_l-y[0]] for x,y in orfs_min_coordinates_rc.items()})
log('ORFS located')



############################ Divide data into ORFs and negatives (non-ORFs)

handle_orf_scores = 'Analysis/ORF_scores.csv'
handle_negatives_scores = 'Analysis/negatives_scores.csv'
(orfs_plus_scored, negatives_plus_scored) = scorer('plus', orfs_plus_coordinates, all_scores)
(orfs_min_out_scored, negatives_min_scored) = scorer('min', orfs_min_coordinates, all_scores)
RET_score_writer(orfs_plus_scored, negatives_plus_scored, orfs_min_out_scored, negatives_min_scored, handle_out_orfs, handle_out_negatives)

log('Positives and Negatives found and made')




######################





Positives, Negatives = data_parser_orfs(handle_orf_scores, handle_negatives_scores)
# Do a iterative PRC analysis to find ideal cut-off parameters
best_filterset, Pr, Rc, d, TP, FP, FDR = Presicion_recall_main(Positives, Negatives)

write_filtered_orfs(handle_orf_scores, best_filterset)

log ('Presicion: ' +str(Pr)+ ' Recall: ' +str(Rc)+ ' True positives: ' +str(TP)+ ' True negatives: ' +str(FP)+ ' False discovery rate: ' +str(FDR))
log('Precicion-Recall analysis finished!')
log('Filtered ORFs written to file!')

# Make histograms for all 5 parameters with drawn cut-off used





histogram_maker(Positives, Negatives, best_filterset)
log('2D histograms made')


#####################

#  Make META plot for annotated genes

META_annotated(Positives ,features_plus, features_min, 'Analysis/Metaplot.csv')
log('Retapamulin METAPLOT made')