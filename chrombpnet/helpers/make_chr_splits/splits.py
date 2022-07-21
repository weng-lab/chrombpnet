import json
import argparse
import os
import random
import pandas as pd
import numpy as np
from chrombpnet.training.utils import data_utils 

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-b", "--bigwig", type=str, required=True, help="Bigwig of tn5 insertions")
    parser.add_argument("-p", "--peaks", type=str, required=True, help="Peaks")
    parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="Non-Peaks")
    parser.add_argument("-il", "--inputlen", type=int, default=2114, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, default=1000, help="Prediction output length")
    parser.add_argument("-j", "--max_jitter", type=int, default=0, help="Maximum jitter applied on either side of region")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Path to store the fold information")
    args = parser.parse_args()

    peak_regions = pd.read_csv(args.peaks, sep='\t', names=NARROWPEAK_SCHEMA)
    nonpeak_regions = pd.read_csv(args.nonpeaks, sep='\t', names=NARROWPEAK_SCHEMA)

    print("Loading Data")

    peak_seqs, peak_cts, peak_coords, nonpeak_seqs, nonpeak_cts, nonpeak_coords, = data_utils.load_data(peak_regions, nonpeak_regions,
                                                                                                        args.genome, args.bigwig,
                                                                                                        args.inputlen, args.outputlen,
                                                                                                        args.max_jitter)
    peak_cts = np.sum(peak_cts, axis=1)
    nonpeak_cts = np.sum(nonpeak_cts, axis=1)

    print(peak_cts.shape)
    print(peak_coords.shape)

    print(nonpeak_cts.shape)
    print(nonpeak_coords.shape)
    
    peak_cts = peak_cts.tolist()
    peak_chroms = peak_coords[:,0].tolist()
    peak_pos = peak_coords[:,1].tolist()

    nonpeak_cts = nonpeak_cts.tolist()
    nonpeak_chroms = nonpeak_coords[:,0].tolist()
    nonpeak_pos = nonpeak_coords[:,1].tolist()

    all_cts = peak_cts + nonpeak_cts
    all_chroms = peak_chroms + nonpeak_chroms
    all_pos = peak_pos + nonpeak_pos

    print("Creating Splits")

    all_df = pd.DataFrame({'chr': all_chroms, 'pos': all_pos, 'cts': all_cts})
    all_df.sort_values(by='cts', inplace=True)

    all_dict = {'fold0': [], 'fold1': [], 'fold2': [], 'fold3': [], 'fold4': []}

    for index,row in all_df.iterrows():
        if index % 10000 == 0:
            print(index)
        test_or_valid = random.choice(['valid', 'test'])
        test_or_valid_fold = random.choice(range(5))
        for fold in range(5):
            if fold != test_or_valid_fold:
                all_dict['fold' + str(fold)].append('train')
            else:
                all_dict['fold' + str(fold)].append(test_or_valid)

    for fold in range(5):
        all_df['fold' + str(fold)] = all_dict['fold' + str(fold)]

    print("Saving Splits")

    all_df.to_csv(args.output_prefix + '.splits.tsv', sep='\t')

if __name__=="__main__":
    main()
    
