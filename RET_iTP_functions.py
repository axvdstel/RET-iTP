#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 09:56:32 2021

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
from RET_iTP_settings import *
import pandas as pd
import matplotlib.pyplot as plt

def genome_l():
        for genome in SeqIO.parse(genome_file, 'fasta'):
            return len(genome.seq)
genome_l=genome_l()

##################General functions
def WIG_writer(out_path, data, genome):
    genome_l = len(genome.seq)
    with open(out_path , 'w') as f:
            f.write('variableStep chrom=%s' %(genome.name) + '\n')
            for x in range(genome_l):
                if data[x] != 0:
                    f.write(str(x+1) + ' ' + str(data[x]) + '\n')

def read_WIG(handle):
    with open(handle, 'r') as f:
            data = {int(x.split(' ')[0])-1:float(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]}
    return data

def revcom(sequence):
    RC = {'G':'C', 'C':'G', 'A':'T', 'T':'A'}
    rc = "".join(str(RC.get(x,'N')) for x in sequence[::-1])
    return rc



#######

def rev_trimmer(sample_path):
    sample_name = sample_path.split('/')[-1]
    command = '%s -Xmx4g in=RAW_reads/%s outm=Trimmed/%s_trimmed.fq out=Trimmed/%s_trimmed_failed.fq literal=%s ktrim=l k=15 mink=11 hdist=2 hdist2=1 mm=f rcomp=f restrictleft=45 ordered minlen=115' %(BBDUK_path, sample_name,sample_name, sample_name, adapter_sequence)
    out = subprocess.run(command, capture_output=True, shell=True, encoding="utf8")
    with open ('log_adapter_trim_%s' %(sample_name) 'w') as f:
        f.write(out.stdout)
    return

def minimapper(sample_path):
    sample_name = sample_path.split('/')[-1]
    command = '%s -x sr %s Trimmed/%s_trimmed.fq > Mapped/%s.paf' %(minimap2_path, genome_file, sample_name, sample_name)
    out = subprocess.run(command, capture_output=True, shell=True, encoding="utf8")
    with open ('log_minimapper_%s' %(sample_name) 'w') as f:
        f.write(out.stdout)
    return



##########

def find_Psites(sample_path):
    if sample_path == library_path:
        return
    sample_name = sample_path.split('/')[-1]
    handle = 'Mapped/%s.paf' %(sample_name)
    counts_plus = {x:0 for x in range(genome_l)}
    counts_min = {x:0 for x in range(genome_l)}
    with open(handle, 'r') as f:
        for line in f:
            strand = line.split('\t')[4]
            left_coor_query = int(line.split('\t')[2])
            left_coor_genome = int(line.split('\t')[7])
            right_coor_genome = int(line.split('\t')[8])
            if strand == '+':
                RET_FP = (left_coor_genome - left_coor_query + adapter_extra + ribosome_protect)%genome_l
                counts_min[RET_FP] += 1
            elif strand == '-':
                RET_FP = (right_coor_genome + left_coor_query - adapter_extra - ribosome_protect - 1)%genome_l
                counts_plus[RET_FP] += 1

    WIG_writer('Analysis/RET-iTP_reads_%s_plus.wig' %(sample_name), counts_plus)
    WIG_writer('Analysis/RET_iTP_reads_%s_min.wig' %(sample_name), counts_min)
    return

def WIG_writer_background():
    library_name = library_path.split('/')[-1]
    with open('Mapped/%s' %(library_name), 'r') as f:
            content = []
            for line in f:
                line_split = line.split('\t')
                content.append([line_split[0], line_split[4], int(line_split[7]), int(line_split[8])])

    counts_plus = {x:0 for x in range(genome_l)}
    counts_min = {x:0 for x in range(genome_l)}
    already_done = 0
    length = len(content)
    for i,line in enumerate(content[:-1]):
        if already_done == 1:
            already_done = 0
            continue
        strand = line[1]
        read_name = line[0]
        read_name_next = content[i+1][0]

        if read_name == read_name_next:
            coor = [line[2], line[3], content[i+1][2], content[i+1][3]]
            left_coor_genome = min(coor)
            right_coor_genome = max(coor)
            already_done = 1

        else:
            left_coor_genome = line[2]
            right_coor_genome = line[3]
        if right_coor_genome - left_coor_genome > 300 :
            continue
        if strand == '+':
            for x in range(left_coor_genome,right_coor_genome):
                counts_plus[x] += 1
        elif strand == '-':
            for x in range(left_coor_genome,right_coor_genome):
                counts_min[x] += 1
        if i%100000 == 0:
            print(i / length)

    WIG_writer('Analysis/background_plus.wig', counts_plus)
    WIG_writer('Analysis/background_min.wig', counts_min)
    return

##############
def background_parser():
    handle_plus = 'Analysis/background_plus.wig'
    handle_min = 'Analysis/background_min.wig'
    with open(handle_plus, 'r') as f:
        bg_plus = np.array([int(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])
    with open(handle_min, 'r') as f:
        bg_min = np.array([int(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])
    return bg_plus, bg_min

def smoother(data):
    data_new=[]
    for i in range(len(data)):
        if i < 150:
            left=i
            right=150
        if i > (len(data)-150):
            left=150
            right = len(data)-i
        else:
            left, right =150,150
        if len(data[i-left:i+right]) > 0 :
            average = sum(data[i-left:i+right]) / len(data[i-left:i+right])
        else:
            average = data[i]
        data_new.append(average)
    return data_new

def bg_lim(data_new):
    data_new = np.array(data_new)
    data_new_mean = np.mean(data_new)
    data_new_range = (data_new_mean -(2 * np.std(data_new)), data_new_mean +(2 * np.std(data_new)))
    data_new = [x if x > data_new_range[0] else data_new_range[0] for x in data_new]
    data_new = [x if x < data_new_range[1] else data_new_range[1] for x in data_new]
    return data_new
############
def frequency_reads(data): #Return data as frequency per 100.000.000 reads
    data = np.array([data.get(x,0) for x in range(genome_l)])
    sum_reads = sum(data)
    norm_data = data / sum_reads *100000000
    return norm_data

##############
def orf_finder(genome):
    seq = genome.seq
    min_orf = 3
    max_orf = 10000
    orfs = defaultdict()
    starts = [codon.start() for codon in re.finditer('ATG|GTG|TTG', seq)]

    for x in starts:
        stops = [codon.start()+x for codon in re.finditer('TAG|TGA|TAA', seq[x:x+max_orf])]

        for y in stops:
            if (y-x) > 0 and (y-x) % 3 == 0:
                if (y-x) < min_orf:
                    break
                if min_orf <= (y-x) <= max_orf:
                    orfs[x]=[y+2]
                    break
                elif (y-x) >max_orf:
                    break
    return orfs


#################

def RET_iTP_scorer():
    #PLUS STRAND
    handle_plus1 = 'Analysis/RET-iTP_scores_%s_plus.wig' %(RET_1_path.split('/')[-1])
    handle_min1 = 'Analysis/RET-iTP_scores_%s_min.wig' %(RET_1_path.split('/')[-1])
    handle_plus_RETmin1 = 'Analysis/RET-iTP_scores_%s_plus.wig' %(control_1_path.split('/')[-1])
    handle_min_RETmin1 = 'Analysis/RET-iTP_scores_%s_plus.wig' %(control_1_path.split('/')[-1])

    handle_plus2 = 'Analysis/RET-iTP_scores_%s_plus.wig' %(RET_2_path.split('/')[-1])
    handle_min2 = 'Analysis/RET-iTP_scores_%s_min.wig' %(RET_2_path.split('/')[-1])
    handle_plus_RETmin2 = 'Analysis/RET-iTP_scores_%s_plus.wig' %(control_2_path.split('/')[-1])
    handle_min_RETmin2 = 'Analysis/RET-iTP_scores_%s_plus.wig' %(control_2_path.split('/')[-1])

    with open(handle_plus1, 'r') as f:
        data_plus1 = np.array([1+float(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])
    with open(handle_min1, 'r') as f:
        data_min1 = np.array([1+float(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])
    with open(handle_plus_RETmin1, 'r') as f:
        data_plus_RETmin1 = np.array([1+float(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])
    with open(handle_min_RETmin1, 'r') as f:
        data_min_RETmin1 = np.array([1+float(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])

    with open(handle_plus2, 'r') as f:
        data_plus2 = np.array([1+float(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])
    with open(handle_min2, 'r') as f:
        data_min2 = np.array([1+float(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])
    with open(handle_plus_RETmin2, 'r') as f:
        data_plus_RETmin2 = np.array([1+float(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])
    with open(handle_min_RETmin2, 'r') as f:
        data_min_RETmin2 = np.array([1+float(x.split(' ')[1]) for x in f.read().split('\n')[1:-1]])

    AverageRETplus = (data_plus1 + data_plus2 ) / 2
    errorRETplus = np.absolute(data_plus1 - data_plus2)
    AverageRETmin = (data_plus_RETmin1 + data_plus_RETmin2) / 2
    errorRETmin = np.absolute(data_plus_RETmin1 - data_plus_RETmin2)
    enrichment = AverageRETplus / AverageRETmin
    enrichment_error = np.sqrt(np.square(errorRETmin/AverageRETmin)+np.square(errorRETplus/AverageRETplus)) * enrichment
    error_frac = errorRETplus / (2*AverageRETplus)

    background_enrichment=[]
    for pos in range(genome_l-1):
        background = [AverageRETplus[x] for x in range(pos-5, (pos+5)%genome_l)]
        try:
            background.remove(AverageRETplus[pos])
            background.remove(max(background))
        except:
            print(pos)
        background_enrichment.append(AverageRETplus[pos] / (sum(background)/8) )
    background_enrichment = np.array(background_enrichment)

    ID = np.array(["%s-plus" %(x) for x in range(genome_l-1)])
    outlist_plus = pd.DataFrame({'AverageRETplus':AverageRETplus, 'errorRETplus':errorRETplus, 'AverageRETmin':AverageRETmin, 'errorRETmin':errorRETmin, 'enrichment':enrichment, 'enrichment_error':enrichment_error, 'error_frac':error_frac,'background_enrichment':background_enrichment },  index=ID)

    AverageRETplus = (data_min1 + data_min2 ) / 2
    errorRETplus = np.absolute(data_min1 - data_min2)
    AverageRETmin = (data_min_RETmin1 + data_min_RETmin2) / 2
    errorRETmin = np.absolute(data_min_RETmin1 - data_min_RETmin2)
    enrichment = AverageRETplus / AverageRETmin
    enrichment_error = np.sqrt(np.square(errorRETmin/AverageRETmin)+np.square(errorRETplus/AverageRETplus)) * enrichment
    error_frac = errorRETplus / (2*AverageRETplus)

    background_enrichment=[]
    for pos in range(genome_l-1):
        background = [AverageRETplus[x] for x in range(pos-5, (pos+5)%genome_l)]
        try:
            background.remove(AverageRETplus[pos])
            background.remove(max(background))
        except:
            print(pos)
        background_enrichment.append(AverageRETplus[pos] / (sum(background)/8) )
    background_enrichment = np.array(background_enrichment)

    ID = np.array(["%s-min" %(x) for x in range(genome_l-1)])
    outlist_min = pd.DataFrame({'AverageRETplus':AverageRETplus, 'errorRETplus':errorRETplus, 'AverageRETmin':AverageRETmin, 'errorRETmin':errorRETmin, 'enrichment':enrichment, 'enrichment_error':enrichment_error, 'error_frac':error_frac,'background_enrichment':background_enrichment },  index=ID)

    out_list = outlist_plus.append(outlist_min)
    out_list.to_csv(handle_out)
    return out_list

def scorer(strand, orfs, freq_reads):

    if strand == 'plus':
        orf_starts_plus = set([x for x in orfs.keys()])
        negatives_plus = defaultdict(lambda:1, {x:[] for x in range(genome_l)})
        for pos in range(genome_l-1):
            try:
                no_go = set([pos+x for x in range(-4,3)])
                if pos -2 in orf_starts_plus:
                    orfs[pos-2].extend([x for x in freq_reads.loc['%s-plus' %(pos)]])

                elif len(no_go.intersection(orf_starts_plus)) > 0:
                    continue
                else:
                    if True in (x in genome[pos-3:pos+9] for x in ['TAG', 'TGA', 'TAA']):
                        continue
                    else:
                        negatives_plus[pos].extend([x for x in freq_reads.loc['%s-plus' %(pos)]])
            except:
                print(pos, strand)
                continue

        return orfs_plus, negatives_plus
    if strand == 'min':
        orf_starts_min = set([x for x in orfs.keys()])
        negatives_min = defaultdict(lambda:1, {x:[] for x in range(genome_l)})
        for pos in range(genome_l-1):
            try:
                no_go = set([pos+x for x in range(-1,6)])
                if pos +3 in orf_starts_min:
                    orfs[pos+3].extend([x for x in freq_reads.loc['%s-min' %(pos)]])
                elif len(no_go.intersection(orf_starts_min)) > 0:
                    continue
                else:
                    if True in (x in revcom(genome[pos-9:pos+3]) for x in ['TAG', 'TGA', 'TAA']):
                        continue
                    else:
                        negatives_min[pos].extend([x for x in freq_reads.loc['%s-min' %(pos)]])
            except:
                 print(pos, strand)
                 continue
        return orfs_min, negatives_min

def RET_score_writer(orfs_plus, negatives_plus, orfs_min, negatives_min, handle_out_orfs, handle_out_negatives):
    with open(handle_out_orfs, 'w') as f:
        f.write('strand,coor_left,coor_right,annotation,ORF,translation,RET_score,RET_error,RETmin_score,RET_min_error,enrichment_score,enrichment_error,error_frac,background_enrichment' + '\n')
        for k, v in orfs_plus.items():
            strand  = 'plus'
            coor_left = k
            coor_right, RET_score, RET_error, RETmin_score, RET_min_error, enrichment_score, enrichment_error, error_frac, background_enrichment = v
            ORF = genome[coor_left:coor_right+1]
            translation = str(Seq(ORF, generic_dna).translate(table=11))
            try:
                annotation = features_plus[(coor_left+1, coor_right+1)][8].split(';')[5][-4:]
            except:
                annotation = ''
            f.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s \n' %(strand, coor_left,coor_right,annotation,ORF, translation,RET_score, RET_error, RETmin_score, RET_min_error, enrichment_score, enrichment_error, error_frac, background_enrichment))
        for k, v in orfs_min.items():
            strand  = 'min'
            coor_right = k
            coor_left, RET_score, RET_error, RETmin_score, RET_min_error, enrichment_score, enrichment_error = v
            ORF = revcom(genome[coor_left-1:coor_right])
            translation = str(Seq(ORF, generic_dna).translate(table=11))
            try:
                annotation = features_min[(coor_left, coor_right)][8].split(';')[5][-4:]
            except:
                annotation = ''

            f.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s \n' %(strand, coor_left,coor_right,annotation,ORF, translation,RET_score, RET_error, RETmin_score, RET_min_error, enrichment_score, enrichment_error, error_frac,background_enrichment ))

        with open(handle_out_negatives, 'w') as f:
        f.write('strand,coor,RET_score,RET_error,RETmin_score,RET_min_error,enrichment_score,enrichment_error,error_frac,background_enrichment' + '\n')

        for k, content in negatives_plus.items():
            strand  = 'plus'
            coor = k
            try:
                RET_score, RET_error, RETmin_score, RET_min_error, enrichment_score, enrichment_error, error_frac, background_enrichment = content
            except:
                continue
            f.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s \n' %(strand, coor,RET_score, RET_error, RETmin_score, RET_min_error, enrichment_score, enrichment_error, error_frac, background_enrichment))
        for k, content in negatives_min.items():
            strand  = 'min'
            coor = k
            try:
                RET_score, RET_error, RETmin_score, RET_min_error, enrichment_score, enrichment_error, error_frac, background_enrichment = content
            except:
                continue
            f.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s \n' %(strand, coor,RET_score, RET_error, RETmin_score, RET_min_error, enrichment_score, enrichment_error, error_frac, background_enrichment))
    return



#############


def META_annotated(freq_reads,features_plus, features_min, handle_out):

    meta_RETplus={x:[] for x in range(40)}
    meta_RETmin={x:[] for x in range(40)}

    features_plus_starts = set([x[0]+1 for x in features_plus.keys()])
    features_min_starts =set([x[1] for x in features_min.keys()])

    for gene in features_plus_starts:
        coor = gene+1 #Get FP peak coordinate (second ucleotide)
        try:
            gene_score_meta = [freq_reads.loc['%s-plus' %(x)].AverageRETplus for x in range(coor-20,coor+20)]
        except:
            print('ERROR, META_PLOT:')
            print(coor, freq_reads.loc['%s-plus' %(coor)])
            continue
        gene_score_RETmin_meta = [freq_reads.loc['%s-plus' %(x)].AverageRETmin for x in range(coor-20,coor+20)]
        for x in range(40):
            meta_RETplus[x].append(gene_score_meta[x])
            meta_RETmin[x].append(gene_score_RETmin_meta[x])
    for gene in features_min_starts:
        coor = gene-2 #Get FP peak coordinate
        try:
            gene_score_meta = [freq_reads.loc['%s-min' %(x)].AverageRETplus for x in range(coor-20,coor+20)][::-1]
        except:
            print('ERROR, META_PLOT:')
            print(coor, freq_reads['%s-min' %(coor)])
            continue
        try:
            gene_score_RETmin_meta = [freq_reads.loc['%s-min' %(x)].AverageRETmin for x in range(coor-1-20,coor-1+20)][::-1]
        except:
            print('ERROR, META_PLOT:')
            print(coor, freq_reads['%s-min' %(coor)])
            continue
        for x in range(40):
            meta_RETplus[x].append(gene_score_meta[x])
            meta_RETmin[x].append(gene_score_RETmin_meta[x])

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    ax.set_title('Metaplot RET-iTP around annotated startcodons')
    ax.set_ylabel('Average score')
    ax.set_xlabel('coordinate')
    plt.scatter(range(-20,20),[np.mean(meta_RETplus[x]) for x in meta_RETplus], '-o' )
    plt.scatter(range(-20,20),[np.mean(meta_RETmin[x]) for x in meta_RETplus], '-o' )
    plt.legend()
    plt.saveimage('Analysis/metaplot.png', dpi=900)
    plt.close()

    with open(handle_out, 'w') as f:
        f.write('coor;RETplus_mean;RETplus_std;RETmin_mean; RETmin_std' + '\n')
        for x in range(40):
            f.write(str(x) + ';' + str(np.mean(meta_RETplus[x])) + ';'+ str(np.std(meta_RETplus[x])) + ';' + str(np.mean(meta_RETmin[x]))+';'+ str(np.std(meta_RETmin[x])) + '\n')
    return

##########################



def calc_Pr_Rc(positives, negatives, filterset):
    import math
    min_RETplus, min_RETmin , min_enrich, min_bg_enrich, max_error_frac = filterset
    P = positives.shape[0]
    orfs_f = positives[(positives.AverageRETplus>min_RETplus)&(positives.AverageRETmin>min_RETmin)&(positives.enrichment>min_enrich)&(positives.background_enrichment>min_bg_enrich)&(positives.error_frac<max_error_frac)]
    TP = orfs_f.shape[0]
    negatives_f = negatives[(negatives.AverageRETplus>min_RETplus)&(negatives.AverageRETmin>min_RETmin)&(negatives.enrichment>min_enrich)&(negatives.background_enrichment>min_bg_enrich)&(negatives.error_frac<max_error_frac)]
    FP = negatives_f.shape[0]

    try:
        Pr = TP/(TP+FP)
    except:
            Pr = 1
    Rc = TP/P
    try:
        d = np.sqrt(math.pow( (1-Pr), 2) + math.pow((1-Rc), 2) )
    except:
        d=1

    # Filter ORFS based on ideal cut-offs and write to file
    if write:
        orfs_f.to_csv(write)
    return Pr, Rc, d, TP, FP



def Presicion_recall_main(Positives, Negatives):
    from operator import itemgetter

    #initial filters:
    RET_plus_mean_Pos = Positives['AverageRETplus'].mean()
    RET_plus_mean_Neg = Negatives['AverageRETplus'].mean()
    min_RETplus = (RET_plus_mean_Pos + RET_plus_mean_Neg) / 2
    RET_min_mean_Pos = Positives['AverageRETmin'].mean()
    RET_min_mean_Neg = Negatives['AverageRETmin'].mean()
    min_RETmin = (RET_min_mean_Pos + RET_min_mean_Neg) / 2
    min_enrich = 1
    min_bg_enrich = 1
    max_error_frac = 1
    RETplus_step =  (RET_plus_mean_Pos - min_RETplus) / 5
    RETmin_step = (RET_min_mean_Pos - min_RETmin) / 5
    enrich_step = 0.5
    bg_enrich_step = 0.5
    error_step = 0.05
    filters = [min_RETplus, min_RETmin , min_enrich, min_bg_enrich, max_error_frac]
    steps = [RETplus_step, RETmin_step, enrich_step, bg_enrich_step, error_step]

    # Do a iterative PRC analysis to find ideal cut-off parameters
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    ax.plot([0,1],[1,0])
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_title('Precision-Recall curve')
    ax.set_ylabel('Precision ( TP/(TP+FP) )')
    ax.set_xlabel('Recall ( TP/P )')
    for _ in range(10):
        old_filters, old_steps = list(filters), list(steps)
        if _ > 0:
            for filter_char in range(5):
                steps[filter_char] = steps[filter_char] * 0.8


        for filter_char in range(5):

            filters_char = [filters[filter_char]+(steps[filter_char]*x) for x in range(-5,5,1)]
            filter_sets =[]
            for step in filters_char:
                temp_filterset = list(filters)
                temp_filterset[filter_char] = step
                filter_sets.append(temp_filterset)

            Pr_Rc_data = []
            for filterset in filter_sets:
                Pr, Rc, d, TP, FP = calc_Pr_Rc(Positives, Negatives, filterset)
                Pr_Rc_data.append((Pr, Rc, d, filterset))
                ax.scatter(Rc, Pr, marker='.', s=4,  c='black', lw=0)
            best_filterset = sorted(Pr_Rc_data, key=itemgetter(2))[0][3]
            filters[filter_char] = best_filterset[filter_char]
            best_data_point = sorted(Pr_Rc_data, key=itemgetter(2))[0]
    #### ROUND UP best filterset
    best_data_point = sorted(Pr_Rc_data, key=itemgetter(2))[0]
    ax.scatter(best_data_point[1], best_data_point[0], marker='*', c='red')

    # Filter ORFS based on ideal cut-offs and write to file

    Pr, Rc, d, TP, FP = calc_Pr_Rc(Positives, Negatives, best_filterset)
    FDR = FP/(FP+TP)
    plt.savefig('Analysis/Presicion-recall_plot.png' , dpi=900)
    plt.close()
    return best_filterset, Pr, Rc, d, TP, FP, FDR

def data_parser_orfs(handle_P, handle_N):
    orfs = pd.read_csv(handle_P)
    Positives = orfs.query('annotation != " "')
    Negatives = pd.read_csv(handle_N)
    return Positives, Negatives

def write_filtered_orfs(handle_orf_scores, filterset):
    filtered_ORF_handle = 'Analysis/filtered_ORFS.csv'
    orfs = pd.read_csv(handle_orf_scores)
    min_RETplus, min_RETmin , min_enrich, min_bg_enrich, max_error_frac = filterset
    orfs_f = orfs[(orfs.AverageRETplus>min_RETplus)&(orfs.AverageRETmin>min_RETmin)&(orfs.enrichment>min_enrich)&(orfs.background_enrichment>min_bg_enrich)&(orfs.error_frac<max_error_frac)]
    orfs_f.to_csv(filtered_ORF_handle)
    return

##############
def histogram_maker(Positives, Negatives,best_filterset):


    bins_list=[0,0.005524272, 0.0078125,	0.011048543,	0.015625,	0.022097087,	0.03125,	0.044194174,	0.0625,	0.088388348,	0.125,	0.176776695,	0.25,	0.353553391,	0.5,	0.707106781,	1,	1.414213562,	2,	2.828427125,	4,	5.656854249, 8, 11.3137085,	16,	22.627417,	32,45.254834,64,90.50966799,128,181.019336,256,362.038672,512,724.0773439,1024,1448.154688,2048,2896.309376,4096,5792.618751,8192]

    for i,filter_char_name in enumerate(['AverageRETplus', 'AverageRETmin', 'enrichment', 'background_enrichment', 'error_frac']):
        RET_reads_histo_P = np.histogram(Positives[filter_char_name], bins=bins_list)
        RET_reads_histo_N = np.histogram(Negatives[filter_char_name], bins=bins_list)

        fig, ax = plt.subplots()
        ax.set_xscale('log', basex=2)
        plt.plot(RET_reads_histo_P[1][1:], RET_reads_histo_P[0]/sum(RET_reads_histo_P[0]), label='Positives; annotated')
        plt.plot(RET_reads_histo_P[1][1:], RET_reads_histo_N[0]/sum(RET_reads_histo_N[0]), label = 'Negatives; non-starts')
        plt.axvline(x=best_filterset[i], color='black', dashes=(1,1))
        plt.xlabel(filter_char_name)
        plt.ylabel('normalised occurences')
        plt.legend()

        handle_out = 'Analysis/histogram_%s' %(filter_char_name)
        plt.savefig(handle_out, dpi=900)
        plt.close()
    return

