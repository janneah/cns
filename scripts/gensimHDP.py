import pandas as pd
from gensim.models import HdpModel
from gensim import corpora
import argparse

parser = argparse.ArgumentParser(description='HDP Analysis')
parser.add_argument('featurefile', help='Input file of discretized features')
parser.add_argument('outputfile', help='Name of the outputfile')

args = parser.parse_args()

file = pd.read_table(args.featurefile, dtype=str)
df = file.drop(['Sample', 'Chr'], axis = 1)
listedDf = df.values.tolist()

dirichlet_dict = corpora.Dictionary(listedDf)
bow_corpus = [dirichlet_dict.doc2bow(text) for text in listedDf]
hdp_model = HdpModel(
    corpus=bow_corpus, 
    id2word=dirichlet_dict,
    random_state=42
    )

hdp_model.save(args.outputfile)