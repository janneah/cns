import argparse, glob
import feature_func as features
import pandas as pd

parser = argparse.ArgumentParser(description='Create feature files')
parser.add_argument('samplefiles', nargs='+', help='Input file of samples')
parser.add_argument('centroinfo', help='File containing centromere information')
parser.add_argument('ascat', help='File containing ploidy information')
parser.add_argument('gccontent', help='File containing GC content')
parser.add_argument('outputfile', help='Path and name of the output')

args = parser.parse_args()

def fullDF(ascat, centroinfo, gccontent, output):
    feature_df = features.makefeatfile(ascat, centroinfo, gccontent)

    # feature_df = features.discretize(feature_df)
    feature_df.to_csv(output, header=True, index=False, sep='\t')

    return feature_df

fullDF(args.samplefiles, args.centroinfo, args.ascat, args.gccontent, args.outputfile)

 
# sample = glob.glob('/home/janneae/TCGA/DerivedData/PanCancer/TCGA_ASCAT_RAW_PVL/ASCAT_TCGA/TCGA-02-*.segments.raw.txt', recursive=True)
# # sample = glob.glob('/home/janneae/TCGA/DerivedData/PanCancer/TCGA_ASCAT_RAW_PVL/ASCAT_TCGA/TCGA-02-0001.segments.raw.txt', recursive=True)
# centromereinfo = '/home/janneae/cns/data/chrominfo.snp6.txt'
# ascat = '/home/janneae/cns/data/filteredAscat.txt'
# gc = '/home/janneae/cns/data/gc.content.txt'
# output = '/home/janneae/cns/steps/test.txt'

# features.raw(sample, ascat)

# full = fullDF(sample, centromereinfo, ascat, gc, output)

