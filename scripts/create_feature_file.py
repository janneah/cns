import argparse
import feature_func as features

parser = argparse.ArgumentParser(description='Create feature files')
parser.add_argument('ascat', help='File containing ploidy information')
parser.add_argument('centroinfo', help='File containing centromere information')
parser.add_argument('gccontent', help='File containing GC content')
parser.add_argument('repeats', help='File containing repeat information')
parser.add_argument('outputfile', help='Path and name of the output')


args = parser.parse_args()

def fullDF(ascat, centroinfo, gccontent, repeats, output):
    feature_df = features.makefeatfile(ascat, centroinfo, gccontent, repeats)
    nbins = [4, 5, 6, 7, 8, 9, 10]
    for i in nbins:
        feature_df = features.discretize(feature_df, i)
        feature_df.to_csv(f'{output}_{i}.features', header=True, index=False, sep='\t')

    return feature_df

fullDF(args.ascat, args.centroinfo, args.gccontent, args.repeats, args.outputfile)