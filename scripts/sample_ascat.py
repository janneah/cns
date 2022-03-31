import pandas as pd
import numpy as np
import os, argparse

parser = argparse.ArgumentParser(description='Take a subset of x samples from each cancertype')
parser.add_argument('ascat', help='Input file of samples')
parser.add_argument('nsamples')
parser.add_argument('output')

args = parser.parse_args()


def sample(df, nsamples):
    cancertypes = np.unique(df['cancer_type'])
    samples = df.groupby(['ID'])
    
    sampled_df = []
    for i in range(len(cancertypes)):
        subset = df[df.cancer_type == cancertypes[i]].sample(n=nsamples, random_state=17).values.tolist()
        for j in subset:
            sampled_df.append(j)
    
    subset_df = pd.DataFrame(sampled_df, columns=df.columns)
    final_df = df[df['ID'].isin(subset_df['ID'])]
    
    return final_df

df = pd.read_table(args.ascat, dtype=str)
new_df = sample(df, args.nsamples)
new_df.to_csv(args.output, header=True, index=False, sep='\t')