import json
import argparse
import os
import random
from chrombpnet.training.utils import data_utils 

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-b", "--bigwig", type=str, required=True, help="Bigwig of tn5 insertions")
    parser.add_argument("-p", "--peaks", type=str, required=True, help="Peaks")
    parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="Non-Peaks")
    parser.add_argument("-il", "--inputlen", type=int, default=2114, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, default=1000, help="Prediction output length")
    parser.add_argument("-j", "--max_jitter", type=int, default=0, help="Maximum jitter applied on either side of region")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Path to store the fold information")
    args = parser.parse_args()

    peak_regions = pd.read_csv(args.peaks, header=None, sep='\t')
    nonpeak_regions = pd.read_csv(args.nonpeaks, header=None, sep='\t')

    peak_seqs, peak_cts, peak_coords, nonpeak_seqs, nonpeak_cts, nonpeak_coords, = data_utils.load_data(peak_regions, nonpeak_regions,
                                                                                                        args.genome, args.bigwig,
                                                                                                        args.inputlen, args.outputlen,
                                                                                                        args.max_jitter)

    peak_df = pd.DataFrame({'coords': peak_coords, 'cts': peak_cts, 'fold0': [], 'fold1': [], 'fold2': [], 'fold3': [], 'fold4': []})
    peak_df.sort_values(by='cts', inplace=True)

    for index,row in peak_df.iterrows():
        test_or_valid = random.choice(['valid', 'test'])
        test_or_valid_fold = random.choice(range(5))
        for fold in range(5):
            if fold != test_or_valid_fold:
                peak_df['fold' + str(fold)].append('train')
            else:
                peak_df['fold' + str(fold)].append(test_or_valid)

    
    peak_df.to_csv(args.outdir + '/splits.tsv', sep='\t')

if __name__=="__main__":
    main()
    
