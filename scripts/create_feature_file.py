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

def fullDF(samplefiles, centroinfo, ascat, gccontent, output):
    sampledfs = features.listdfs(samplefiles)
    lst = []
    for i in sampledfs:
        lst.append(features.makefeatfile(i, centroinfo, ascat, gccontent))

    feature_df = pd.concat(lst)
    feature_df.to_csv(output, header=True, index=False, sep='\t')

    return feature_df

fullDF(args.samplefiles, args.centroinfo, args.ascat, args.gccontent, args.outputfile)

 
# sample = glob.glob('/home/janneae/TCGA/DerivedData/PanCancer/TCGA_ASCAT_RAW_PVL/ASCAT_TCGA/TCGA-*.segments.raw.txt', recursive=True)
# centromereinfo = '../data/chrominfo.snp6.txt'
# ascat = '../data/filteredAscat.txt'
# gc = '../data/gc_content.txt'
# output = '../steps/testall_features_mateo10.txt'

# fullDF(sample, centromereinfo, ascat, gc, output)


