import csv
import re
import pandas as pd
import collections
import numpy as np
from scipy import stats
from math import trunc
import threading
import multiprocessing as mp
from Bio import SeqIO
import pysam
import scipy
import math
import os
import subprocess


path = ''
sample = ''

input_folder = (path + sample) # path of the splited bam files
input_files = [f for f in os.listdir(input_folder) if f.endswith('.bam')]

bam_split = pysam.AlignmentFile((input_folder + '/' + 'Barcode_Title' + '.bam'), 'rb')

Read_CB_UMI = open((input_folder + '/process/' + 'Barcode_Title' + '.txt'),'w') 
Read_CB_UMI.write('ReadsName\tchr\tCellBarcode\tUMI\n')

for read in bam_split:
    try:
        ReadsName = read.query_name
        tags = read.get_tags()
        CellBarcode = [tag[1] for tag in tags if tag[0] == 'DB']
        UMI = [tag[1] for tag in tags if tag[0] == 'UR']
        if len(CellBarcode) != 0 and len(UMI) != 0:
            Read_CB_UMI.write(ReadsName + '\t' + read.reference_name + '\t' + CellBarcode[0] + '\t' + UMI[0] + '\n')
    except:
        continue
Read_CB_UMI.close()

mpileup_split = pd.read_table((input_folder + '/' + 'Barcode_Title' + '.mpileup'),\
    names = ['pos', 'ref','depth','alt','quality','?','ReadsName'])
mpileup_split = mpileup_split[mpileup_split['depth'] != '0'].reset_index().drop('?', axis = 1)
mpileup_split.columns = ['chr','pos','ref','depth','alt','quality','ReadsName']

mpileup_split_HaveMut = mpileup_split[mpileup_split['alt'].str.contains('[a-zA-Z]')].reset_index(drop = True)
#mpileup_split_HaveMut = mpileup_split[mpileup_split['chr'] == 'Plant__ZeawithTE_8'].reset_index(drop = True)

def replace_with_pipe(match):
    before, n, after = match.groups()
    n = int(n)
    return before + str(n) + after[:n] + '|' + after[n:]

mpileup_split_HaveMut_reads_unique = pd.DataFrame(columns=['chr', 'pos', 'ref','alt','ReadsName'])

for i in range(len(mpileup_split_HaveMut['chr'])):
    ref_chr = mpileup_split_HaveMut['chr'][i]
    ref = mpileup_split_HaveMut['ref'][i].upper()
    pos = mpileup_split_HaveMut['pos'][i]
    alt = mpileup_split_HaveMut['alt'][i].upper()

    pattern_1 = r'(\D)(\d+)([acgtnACGTN]+)'
    alt = re.sub(pattern_1, replace_with_pipe, alt)

    alt = re.sub('\^.','',alt)
    alt = re.sub('\$','',alt)
    pattern = r'([,.][+-]\d+[ACGTNacgtn]+)|([,.<>*])|([ACGTNacgtn])'
    matches = re.findall(pattern, alt)
    alt_list = [item for tup in matches for item in tup if item != '']
    ReadsName_list = mpileup_split_HaveMut['ReadsName'][i].split(',')
    for j in range(len(ReadsName_list)):
        mpileup_split_HaveMut_reads_unique = mpileup_split_HaveMut_reads_unique.append({'chr':ref_chr,'pos':pos,'ref': ref,\
            'alt': alt_list[j],'ReadsName':ReadsName_list[j]}, ignore_index=True)

Read_CB_UMI = pd.read_table((input_folder + '/process/' + 'Barcode_Title' + '.txt')) 

mpileup_split_HaveMut_reads_unique = pd.merge(mpileup_split_HaveMut_reads_unique,Read_CB_UMI,on = ['ReadsName','chr'])
mpileup_split_HaveMut_reads_unique = mpileup_split_HaveMut_reads_unique.drop('ReadsName', axis = 1)
mpileup_split_HaveMut_reads_unique.to_csv((input_folder + '/process/' + 'Barcode_Title' + '_MutReads.txt'), sep = '\t')

mpileup_split['CellBarcode'] = 'Barcode_Title'
mpileup_split.to_csv((input_folder + '/process/' + 'Barcode_Title' + '_PosReadsDepth.txt'), sep = '\t')