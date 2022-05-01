import pandas as pd

ascat = '~/cns/data/filteredAscatRaw.txt'
repeats = '~/cns/data/repeats.txt'

def no_repeats(ascat, line1):  
    reps = [0] * len(ascat)
    line1 = line1[line1['repFamily']=='L1'].reset_index()
    line1['Chr'] = line1['genoName'].str.split('_').str[0].str.extract('(\d+)', expand=False).dropna().astype(int)

    for i in range(0, len(ascat)):
        count = 0
        for j in range(0, len(line1)):
            if line1['Chr'][j] == ascat['Chr'][i] and line1['genoStart'][j] >= ascat['Start'][i] and line1['genoStart'][j] < ascat['End'][i]:
                count += 1
            elif line1['Chr'][j] == ascat['Chr'][i] and line1['genoStart'][j] < ascat['Start'][i] and line1['genoEnd'][j] > ascat['Start'][i]:
                count += 1
            else:
                continue
        reps[i] = count
    
    return reps

no_repeats(ascat, repeats)