import math
import pandas as pd
import numpy as np

def relabelSC(df):
    df.loc[df['Chr'] == 'X', 'Chr'] = 23
    df.loc[df['Chr'] == 'Y', 'Chr'] = 24
    return df

def getDist2Centromere(df, anno):
    dist_vec = [0] * len(df)
    for i in range(len(df)):
        segstart = df['Start'][i]
        segend = df['End'][i]
        chrom_no = int(df['Chr'][i])
        
        centstart = anno.loc[anno['chrom'] == chrom_no]['centstart'].iloc[0]
        centend = anno.loc[anno['chrom'] == chrom_no]['centend'].iloc[0]
        
        if segend < centstart:
            dist_vec[i] = centstart - segend
        elif segstart > centend:
            dist_vec[i] = segstart - centend
        else:
            dist_vec[i] = 0
            
    return dist_vec

def getSegVal(ascat):
    SegValVec = [0] * len(ascat)
    
    for i in range(len(ascat)):
        segval = (ascat['nAraw'][i] + ascat['nBraw'][i]) / ascat['Ploidy'][i]
    
        if segval > 0: SegValVec[i] = math.log2(segval)
        else: SegValVec[i] = 0

    return np.around(np.array(SegValVec),2)

def getLOH(df):
    LOH = [0] * len(df)
    
    for i in range(len(df)):
        if df['nA'][i] == 0 or df['nB'][i] == 0: LOH[i] = 1
        else: LOH[i] = 0
    
    return LOH

def getCP(df):
    t = df.groupby(['ID', 'Chr'])['Chr'].transform('count') - 1
    return df.groupby(['ID', 'Chr'])['Chr'].transform('count') - 1

def getSizeofDiploidSeg(df):
    size_diploid = [0] * len(df)

    for i in range(len(df)):
        if df['nA'][i] + df['nB'][i] == 2: size_diploid[i] = df['End'][i] - df['Start'][i]
        else: size_diploid[i] = 0

    return size_diploid

def getGCcontent(df, gc_file):
    gc = [0] * len(df)

    for i in range(len(df)):
        chrom = int(df['Chr'][i])
        chr_rows = gc_file.loc[(gc_file['Chr'] == chrom) & (gc_file['Position'] >= df['Start'][i]) & (gc_file['Position'] <= df['End'][i])]
        gc_bases = chr_rows.loc[chr_rows['base'].isin(['G', 'C'])]['X100bp']
        
        if gc_bases.empty:
            gc[i] = 0
        else:
            gc[i] = gc_bases.mean().round(2)
        
    return gc

def getRepeats(df):
    return

def getDist2CNV(df):
    dist2CNV = [0] * len(df)
    dist2CNV[0] = int(df['Start'][1] - df['End'][0])
    dist2CNV[len(df) - 1] = int(df['Start'][len(df) - 1] - df['End'][len(df) - 2])

    for i in range(1, len(df)-1):
        if df['Chr'][i] == df['Chr'][i - 1]:
            prevCNV = df['End'][i - 1]
            dist2prev =  df['Start'][i] - prevCNV
        elif df['Chr'][i] != df['Chr'][i - 1]:
            dist2prev = float('inf')
        
        if df['Chr'][i] == df['Chr'][i + 1]:
            nextCNV = df['Start'][i + 1]
            dist2next = nextCNV - df['End'][i]
        elif df['Chr'][i] != df['Chr'][i + 1]:
            dist2next = float('inf')
        
        dist2CNV[i] = min(dist2prev, dist2next)

    return dist2CNV

def no_repeats(ascat, line1):
    reps = [0] * len(ascat)
    line1 = line1[line1['repFamily']=='L1'].reset_index()
    line1['Chr'] = line1['genoName'].str.split('_').str[0].str.extract('(\d+)', expand=False)
    line1 = line1.dropna(subset=['Chr'])
    line1 = line1.astype({"Chr": int}, errors='raise') 
    
    for i in range(0, len(ascat)):
        count = 0
        line1s = line1[
            (line1.Chr == ascat.Chr[i]) & (line1.genoStart >= ascat.Start[i]) & (line1.genoStart < ascat.End[i]) 
            | 
            (line1.Chr == ascat.Chr[i]) & (line1.genoStart < ascat.Start[i]) & (line1.genoEnd > ascat.Start[i])
            ]
        
        reps[i] = len(line1s.index)
    
    return reps

def makefeatfile(ascat, centromere, gccontent, repeats):
    centromere = pd.read_table(centromere, sep=' ')
    ascat = pd.read_table(ascat, sep='\t')
    gccontent = pd.read_table(gccontent, sep=' ')
    repeats = pd.read_table(repeats, sep='\t')

    feature_df = pd.DataFrame({
              'Sample': ascat['ID'],
                 'Chr': ascat['Chr'],
                  'CN': ascat['cn'], 
             'SegSize': ascat['End'] - ascat['Start'],
           'Dist2Cent': getDist2Centromere(ascat, centromere),
              'SegVal': getSegVal(ascat),
                 'LOH': getLOH(ascat),
          'SizeDipSeg': getSizeofDiploidSeg(ascat),
                'CpCN': getCP(ascat),
           'Dist2nCNV': getDist2CNV(ascat),
              'GCcSeg': getGCcontent(ascat, gccontent),
           'NumRepeats': no_repeats(ascat, repeats)
     
    })

    feature_df = feature_df.groupby(['Sample', 'Chr'], as_index=False, sort=False).mean()
    feature_df = feature_df[feature_df.Chr != 23]
 
    return feature_df

def discretize(inputdf, nbins):
    labels = [i for i in range(1, nbins + 1)]
    df = inputdf

    df['CN'] = pd.cut(
                        x=df['CN'], 
                        bins=[0, 1, 2, 3, 4, 5, df['CN'].max()], 
                        labels=[1, 2, 3, 4, 5, 6],
                        include_lowest = True
                        )
    df['SegSize'] = pd.cut(
                        x=df['SegSize'], 
                        bins=[0, 1e5, 1e6, 3e6, 1e7, 5e7, 1e1000], 
                        labels=[1, 2, 3, 4, 5, 6],
                        include_lowest=True
                        )
    df['Dist2Cent'] = pd.qcut(
                        x=df['Dist2Cent'], 
                        q=nbins,
                        labels=labels
                        )
    df['SegVal'] = pd.qcut(
                        x=df['SegVal'],
                        q=nbins,
                        labels=labels
                        )
    df['LOH'] = df['LOH'].round().astype(int)
    df['SizeDipSeg'] = pd.cut(
                        x=df['SizeDipSeg'],
                        bins=[0, 1, 5e7, 1e8, 3e8],
                        labels=[1, 2, 3, 4],
                        include_lowest=True
                        )
    df['CpCN'] = pd.qcut(
                        x=df['CpCN'],
                        q=nbins,
                        labels=labels
                        )
    df['Dist2nCNV'] = pd.qcut(
                        x=df['Dist2nCNV'], 
                        q=nbins,
                        labels=labels
                        )
    df['GCcSeg'] = pd.qcut(
                        x=df['GCcSeg'], 
                        q=nbins,
                        labels=labels
                        )
    df['NumRepeats'] = pd.qcut(
                        x=df['NumRepeats'], 
                        q=nbins,
                        labels=labels
                        )
    
    df['CN'] = 'CN_' + df['CN'].astype(str)
    df['SegSize'] = 'SegSize_' + df['SegSize'].astype(str)
    df['Dist2Cent'] = 'Dist2Cent_' + df['Dist2Cent'].astype(str)
    df['SegVal'] = 'SegVal_' + df['SegVal'].astype(str)
    df['LOH'] = 'LOH_' + df['LOH'].astype(str)
    df['SizeDipSeg'] = 'SizeDipSeg_' + df['SizeDipSeg'].astype(str)
    df['CpCN'] = 'CpCN_' + df['CpCN'].astype(str)
    df['Dist2nCNV'] = 'Dist2CNV_' + df['Dist2nCNV'].astype(str)
    df['GCcSeg'] = 'GCcSeg_' + df['GCcSeg'].astype(str)
    df['NumRepeats'] = 'NumRepeats_' + df['NumRepeats'].astype(str)

    return df
