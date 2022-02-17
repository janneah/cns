import math, os
import pandas as pd
import numpy as np

def relabelSC(df):
    df.loc[df['chr'] == 'X', 'chr'] = 23
    df.loc[df['chr'] == 'Y', 'chr'] = 24
    return df

def getDist2Centromere(df, anno):
    dist_vec = [0] * len(df)
    for i in range(len(df)):
        segstart = df['startpos'][i]
        segend = df['endpos'][i]
        chrom_no = int(df['chr'][i])
        
        centstart = anno.loc[anno['chrom'] == chrom_no]['centstart'].iloc[0]
        centend = anno.loc[anno['chrom'] == chrom_no]['centend'].iloc[0]
        
        if segend < centstart:
            dist_vec[i] = centstart - segend
        elif segstart > centend:
            dist_vec[i] = segstart - centend
        else:
            dist_vec[i] = 0
            
    return dist_vec

def getSegVal(df, ascat):
    SegValVec = [0] * len(df)
    
    for i in range(len(df)):
        ploidy = ascat.loc[ascat['ID'] == df['sample'][i]]['Ploidy'].iloc[0]
        segval = (df['nAraw'][i] + df['nBraw'][i]) / ploidy

        if segval > 0: SegValVec[i] = math.log2(segval)
        else: SegValVec[i] = 0

    return np.around(np.array(SegValVec),2)

def getLOH(df):
    LOH = [0] * len(df)
    
    for i in range(len(df)):
        if df['nMajor'][i] == 0 or df['nMinor'][i] == 0: LOH[i] = 1
        else: LOH[i] = 0
    
    return LOH

def getCP(df):
    return df.groupby('chr')['chr'].transform('count') - 1

def getSizeofDiploidSeg(df):
    size_diploid = [0] * len(df)

    for i in range(len(df)):
        if df['nMajor'][i] + df['nMinor'][i] == 2: size_diploid[i] = df['endpos'][i] - df['startpos'][i]
        else: size_diploid[i] = 0

    return size_diploid

def getGCcontent(df, gc_file):
    gc = [0] * len(df)

    for i in range(len(df)):
        chrom = int(df['chr'][i])
        chr_rows = gc_file.loc[(gc_file['Chr'] == chrom) & (gc_file['Position'] >= df['startpos'][i]) & (gc_file['Position'] <= df['endpos'][i])]
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
    dist2CNV[0] = int(df['startpos'][1] - df['endpos'][0])
    dist2CNV[len(df) - 1] = int(df['startpos'][len(df) - 1] - df['endpos'][len(df) - 2])

    for i in range(1, len(df)-1):
        if df['chr'][i] == df['chr'][i - 1]:
            prevCNV = df['endpos'][i - 1]
            dist2prev =  df['startpos'][i] - prevCNV
        elif df['chr'][i] != df['chr'][i - 1]:
            dist2prev = float('inf')
        
        if df['chr'][i] == df['chr'][i + 1]:
            nextCNV = df['startpos'][i + 1]
            dist2next = nextCNV - df['endpos'][i]
        elif df['chr'][i] != df['chr'][i + 1]:
            dist2next = float('inf')
        
        dist2CNV[i] = min(dist2prev, dist2next)

    return dist2CNV

def listdfs(files):
    dfs = []
    for file in files:
        if os.path.getsize(file) > 1:
            dfs.append(pd.read_table(file))
    return dfs

def makefeatfile(df, centromere, ascat, gccontent):
    centromere = pd.read_table(centromere, sep=' ')
    ascat = pd.read_table(ascat, sep=' ')
    gccontent = pd.read_table(gccontent, sep=' ')
    df = relabelSC(df)

    feature_df = pd.DataFrame({
              'Sample': df['sample'],
                 'Chr': df['chr'],
          'CopyNumber': df['nMajor'] + df['nMinor'], 
         'SegmentSize': df['endpos'] - df['startpos'],
     'Dist2Centromere': getDist2Centromere(df, centromere),
              'SegVal': getSegVal(df, ascat),
                 'LOH': getLOH(df),
      'SizeDiploidSeg': getSizeofDiploidSeg(df),
       'ChangepointCN': getCP(df),
     'Dist2NearestCNV': getDist2CNV(df),
        'GCcontentSeg': getGCcontent(df, gccontent)
     
    })

    feature_df = feature_df.groupby(['Sample', 'Chr'], as_index=False, sort=False).mean()
 
    return feature_df
