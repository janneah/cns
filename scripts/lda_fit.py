import pandas as pd
from sklearn.decomposition import LatentDirichletAllocation
import argparse

parser = argparse.ArgumentParser(description='LDA Analysis')
parser.add_argument('featurefile', help='Input file of discretized features')
parser.add_argument('ncomponents', type=int, help='An integer specifying the number of components to use')
parser.add_argument('outputfile', help='Name of the outputfile')

args = parser.parse_args()

def lda_fit(file, n_components, output):
    input = pd.read_table(file)

    # Removing columns with sample name and chromosome
    X = input.drop(['Sample', 'Chr'], axis = 1)
    
    # Fitting the model using the features
    model = LatentDirichletAllocation(
            n_components=n_components,
            max_iter=5,                     # Number of epochs
            learning_method="online",       # Options are batch or online (online is faster, uses mini-batch)
            learning_offset=50.0,           # 
            random_state=0
            )

    fit = model.fit_transform(X)
    df = pd.DataFrame(fit)
    df.to_csv(output, header=False, index=False, sep='\t')
    
    return fit

lda_fit(args.featurefile, args.ncomponents, args.outputfile)

