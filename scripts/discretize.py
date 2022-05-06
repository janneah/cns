import argparse
import pandas as pd
import feature_func as features

parser = argparse.ArgumentParser(description='Discretize features')
parser.add_argument('featurefile')
parser.add_argument('outputfile', help='Path and name of the output')

args = parser.parse_args()

nbins = [4, 5, 6, 7, 8]
for i in nbins:
    featuredf = pd.read_table(args.featurefile, sep='\t')
    discretedf = features.discretize(featuredf, i)
    discretedf.to_csv(f'{args.outputfile}_{i}.features', header=True, index=False, sep='\t')