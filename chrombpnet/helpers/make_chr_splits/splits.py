from itertools import count
import json
import argparse
import os
import random
import pandas as pd
import numpy as np
from chrombpnet.training.utils import data_utils 
import pyBigWig

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-b", "--bigwig", type=str, required=True, help="Bigwig of tn5 insertions")
    parser.add_argument("-p", "--peaks", type=str, required=True, help="Peaks")
    parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="Non-Peaks")
    parser.add_argument("-il", "--inputlen", type=int, default=2114, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, default=1000, help="Prediction output length")
    parser.add_argument("-j", "--max_jitter", type=int, default=500, help="Maximum jitter applied on either side of region")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Path to store the fold information")
    args = parser.parse_args()

    peak_regions = pd.read_csv(args.peaks, sep='\t', names=NARROWPEAK_SCHEMA)
    nonpeak_regions = pd.read_csv(args.nonpeaks, sep='\t', names=NARROWPEAK_SCHEMA)

    print("Loading Data")

    #peak_seqs, peak_cts, peak_coords, nonpeak_seqs, nonpeak_cts, nonpeak_coords, = data_utils.load_data(peak_regions, nonpeak_regions,
    #                                                                                                    args.genome, args.bigwig,
    #                                                                                                    args.inputlen, args.outputlen,
    #                                                                                                    args.max_jitter)
    #peak_cts = np.sum(peak_cts, axis=1)
    #nonpeak_cts = np.sum(nonpeak_cts, axis=1)

    #print(peak_cts.shape)
    #print(peak_coords.shape)

    #print(nonpeak_cts.shape)
    #print(nonpeak_coords.shape)
    
    #peak_cts = peak_cts.tolist()
    peak_chroms = peak_regions['chr']
    peak_pos = peak_regions['start'] + peak_regions['summit']

    #nonpeak_cts = nonpeak_cts.tolist()
    nonpeak_chroms = nonpeak_regions['chr']
    nonpeak_pos = nonpeak_regions['start'] + nonpeak_regions['summit']

    #all_cts = peak_cts + nonpeak_cts
    all_chroms = peak_chroms.tolist() + nonpeak_chroms.tolist()
    all_pos = peak_pos.tolist() + nonpeak_pos.tolist()

    print("Creating Splits")

    all_df = pd.DataFrame({'chr': all_chroms, 'pos': all_pos})
    all_df.sort_values(by=['chr', 'pos'], inplace=True)

    group_dict = {}

    cur_chrom = ''
    cur_group = ''
    last_pos = 0
    for index,row in all_df.iterrows():
        if cur_chrom != '':
            if row['chr'] != cur_chrom:
                cur_chrom = row['chr']
                cur_group += 1
                group_dict[cur_group] = [(row['chr'], row['pos'])]
            else:
                if row['pos'] <= int(last_pos) + int(args.inputlen) + int(2 * args.max_jitter):
                    group_dict[cur_group].append((row['chr'], row['pos']))
                else:
                    cur_group += 1
                    group_dict[cur_group] = [(row['chr'], row['pos'])]
        else:
            cur_chrom = row['chr']
            cur_group = 0
            group_dict[cur_group] = [(row['chr'], row['pos'])]
        last_pos = row['pos']

    groups = []
    group_counts = []

    bw = pyBigWig.open(args.bigwig)

    for group in group_dict:
        groups.append(group)
        sum = 0
        for element in group_dict[group]:
            labels = bw.values(element[0], int(element[1] - (args.inputlen // 2)), int(element[1] + (args.inputlen // 2)))
            labels = np.array(labels)
            labels = np.nan_to_num(labels)
            labels = np.sum(labels)
            sum += labels
        group_counts.append(sum)

    group_df = pd.DataFrame({'groups': groups, 'group_counts': group_counts})
    group_df.sort_values(by='group_counts', inplace=True)
    group_fold_dict = {'fold0': [], 'fold1': [], 'fold2': [], 'fold3': [], 'fold4': []}

    count = 0
    valid_used = []

    for index,row in group_df.iterrows():
        if index % 10000 == 0:
            print(index)
        if count % 2 == 0:
            test_or_valid = 'valid'
        else:
            test_or_valid = 'test'
        test_or_valid_fold = random.choice([i for i in range(5) if i not in valid_used])
        for fold in range(5):
            if fold != test_or_valid_fold:
                group_fold_dict['fold' + str(fold)].append('train')
            else:
                group_fold_dict['fold' + str(fold)].append(test_or_valid)
        count += 1
        valid_used.append(test_or_valid_fold)
        if len(valid_used) == 5:
            valid_used = []

    for fold in range(5):
        group_df['fold' + str(fold)] = group_fold_dict['fold' + str(fold)]

    all_dict = {'chr': [], 'pos': [], 'fold0': [], 'fold1': [], 'fold2': [], 'fold3': [], 'fold4': []}

    for index,row in group_df.iterrows():
        for element in group_dict[row['groups']]:
            all_dict['chr'].append(element[0])
            all_dict['pos'].append(element[1])
            all_dict['fold0'].append(row['fold0'])
            all_dict['fold1'].append(row['fold1'])
            all_dict['fold2'].append(row['fold2'])
            all_dict['fold3'].append(row['fold3'])
            all_dict['fold4'].append(row['fold4'])

    splits_df = pd.DataFrame(all_dict)

    print("Saving Splits")

    splits_df.to_csv(args.output_prefix + '.splits.tsv', sep='\t')

if __name__=="__main__":
    main()
    
