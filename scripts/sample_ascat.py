import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Take a subset of x samples from each cancertype')
parser.add_argument('ascat', help='Input file of samples')
parser.add_argument('nsamples')
parser.add_argument('output')

args = parser.parse_args()

def sample(df, nsamples):
    cancertypes = np.unique(df['cancer_type'])
    new_df = df.drop_duplicates(subset='ID')
    sampled_df = []

    for i in range(len(cancertypes)):
        subset = list(new_df[new_df.cancer_type == cancertypes[i]].sample(frac=float(nsamples), replace=False, random_state=17)['ID'])
        sampled_df.extend(subset)
    
    final_df = df[df['ID'].isin(sampled_df)]
    
    return final_df

df = pd.read_table(args.ascat, dtype=str)
new_df = sample(df, args.nsamples)
new_df.to_csv(args.output, header=True, index=False, sep='\t')
