import pandas as pd
from gensim.models import LdaModel
from gensim import corpora
import argparse

parser = argparse.ArgumentParser(description="LDA Analysis using Python's gensim")
parser.add_argument('featurefile', help='Input file of discretized features')
parser.add_argument('ntopics', type=int, help='An integer specifying the maximum number of topics')
parser.add_argument('outputfile', help='Name of the outputfile')

args = parser.parse_args()

file = pd.read_table(args.featurefile, dtype=str)
df = file.drop(['Sample', 'Chr', 'NumRepeats'], axis = 1)
listedDf = df.values.tolist()

dirichlet_dict = corpora.Dictionary(listedDf)
bow_corpus = [dirichlet_dict.doc2bow(text) for text in listedDf]

modelLDA = LdaModel(
        corpus=bow_corpus,
        id2word=dirichlet_dict,
        num_topics=args.ntopics,
        chunksize=len(bow_corpus),
        passes=10,
        alpha='asymmetric',
        eta='symmetric',
        random_state=42)

modelLDA.save(args.outputfile)

