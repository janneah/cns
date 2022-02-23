import pandas as pd
from sklearn.decomposition import LatentDirichletAllocation
import argparse

parser = argparse.ArgumentParser(description='Create feature files')
parser.add_argument('featurefile', help='Input file of discretized features')
parser.add_argument('ncomponents', type=int, help='An integer specifying the number of components to use')

args = parser.parse_args()

def lda_fit(df, n_components):
    # Removing columns with sample name and chromosome
    X = df.drop(['Sample', 'Chr'], axis = 1)
    
    # Fitting the model using the features
    model = LatentDirichletAllocation(
            n_components=n_components,
            max_iter=5,                     # Number of epochs
            learning_method="online",       # Options are batch or online (online is faster, uses mini-batch)
            learning_offset=50.0,           # 
            random_state=0
            )

    fit = model.fit(X)
    
    return fit

file = pd.read_table(args.featurefile)

lda_fit(file, args.ncomponents)

