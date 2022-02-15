import math
import pandas as pd

def relabelSC(col):
    for i in range(0, len(col)):
        if col[i] == 'X': col[i] = 23
        if col[i] == 'Y': col[i] = 24
    return col

def getDist2Centromere(df, anno):
    dist_vec = [0] * len(df)
    for i in range(0, len(df)):
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
    
    for i in range(0, len(df)):
        ploidy = ascat.loc[ascat['ID'] == df['sample'][i]]['Ploidy'].iloc[0]
        segval = (df['nAraw'][i] + df['nBraw'][i]) / ploidy

        if segval > 0: SegValVec[i] = math.log2(segval)
        else: SegValVec[i] = 0

    return SegValVec

def getLOH(df):
    LOH = [0] * len(df)
    
    for i in range(0, len(df)):
        if df['nMajor'][i] == 0 or df['nMinor'][i] == 0: LOH[i] = 1
        else: LOH[i] = 0
    
    return LOH

def getCP(df):
    return df.groupby('chr')['chr'].transform('count') - 1

def getSizeofDiploidSeg(df):
    return

testfile = pd.read_table("~/TCGA/DerivedData/PanCancer/TCGA_ASCAT_RAW_PVL/ASCAT_TCGA/TCGA-02-0001.segments.raw.txt")
# print(getCP(testfile)[1])