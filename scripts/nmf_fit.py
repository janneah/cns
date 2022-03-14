import pandas as pd
from sklearn.decomposition import NMF
import argparse

parser = argparse.ArgumentParser(description='NMF Analysis')
parser.add_argument('featurefile', help='Input file of discretized features')
parser.add_argument('ncomponents', type=int, help='An integer specifying the number of components to use')
parser.add_argument('outputfile', help='Name of the outputfile')

args = parser.parse_args()

def nmf_fit(file, n_components, output):
    input = pd.read_table(file)

    # Removing columns with sample name and chromosome
    X = input.drop(['Sample', 'Chr'], axis = 1)
    
    # Fitting the model using the features
    model = NMF(
        n_components=n_components, 
        beta_loss="kullback-leibler",
        solver="mu",
        max_iter=1000,
        l1_ratio=0.5,
        random_state=0
        )
    
    fit = model.fit_transform(X)
    df = pd.DataFrame(fit)
    df.to_csv(output, header=False, index=False, sep='\t')
    
    return fit

nmf_fit(args.featurefile, args.ncomponents, args.outputfile)
