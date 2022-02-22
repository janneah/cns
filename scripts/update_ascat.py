import pandas as pd
import os, argparse

parser = argparse.ArgumentParser(description='Create feature files')
parser.add_argument('samplefiles', nargs='+', help='Input file of samples')
parser.add_argument('ascat')
parser.add_argument('output')

args = parser.parse_args()

def listdfs(files):
    dfs = []
    for file in files:
        if os.path.getsize(file) > 1:
            dfs.append(pd.read_table(file))
    return dfs

def raw(samples, ascat, output):
    ascat = pd.read_table(ascat, sep=' ')
    sampledfs = listdfs(samples)
    ascat['nAraw'], ascat['nBraw'] = 0, 0
    for i in range(len(ascat)):
        for df in sampledfs:
            if df['sample'][0] == ascat['ID'][i]:
                nAraw = df.loc[(df['startpos'] == ascat['Start'][i]) & (df['endpos'] == ascat['End'][i])]['nAraw'].iloc[0]
                nBraw = df.loc[(df['startpos'] == ascat['Start'][i]) & (df['endpos'] == ascat['End'][i])]['nBraw'].iloc[0]  
                ascat['nAraw'][i] = nAraw
                ascat['nBraw'][i] = nBraw
                break
    ascat.to_csv(output, header=True, index=False, sep='\t')
    return ascat

raw(args.samplefiles, args.ascat, args.output)