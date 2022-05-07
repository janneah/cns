from gwf import Workflow, AnonymousTarget

gwf = Workflow()

"""
Workflow for LDA analysis of ASCAT data
"""
# Various input files
sample = '/home/janneae/TCGA/DerivedData/PanCancer/TCGA_ASCAT_RAW_PVL/ASCAT_TCGA/TCGA*.segments.raw.txt' 
centromereinfo = '../data/chrominfo.snp6.txt'
ascat = '../data/filteredAscat.txt'
gc = '../data/gc.content.txt'
repeats = '../data/repeats.txt'

# Various parameters
nfeat = 10
start = 1
ntopics = 15
nsamples = 0.7

# Intermediary files
updatedascat = '../data/filteredAscatRaw.txt'
sampledascat = f'../steps/sampledascat/sampled_{nsamples}.ascat'
featurefile = f'../steps/featurefiles/discretized_{nfeat}_{nsamples}.features'
discrete_features = f'../steps/featurefiles/discretized_{nfeat}_{nsamples}'

# Validation input and output
remaining30ascat = '../steps/sampledascat/sampled_0.3.ascat'
validation_featurefile = f'../steps/featurefiles/nondiscretized_{nfeat}_0.3.features'

def update_ascat(samplefiles, ascat, updatedascat):
    inputs = [ascat]
    outputs = [updatedascat]
    options = {
        'memory': '8g',
        'walltime': '4-00:00:00',
        'account': 'CancerEvolution'
    }
    spec = f'''
    
    python update_ascat.py {samplefiles} {ascat} {updatedascat}

    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def sample_ascat(ascat, nsamples, output):
    inputs = [ascat]
    outputs = [output]
    options = {
        'memory': '8g',
        'walltime': '02:00:00',
        'account': 'CancerEvolution'
    }
    spec = f'''
    
    python sample_ascat.py {ascat} {nsamples} {output}

    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def create_feature_file(ascat, centromere, gc, repeats, output): 
    inputs = [centromere, ascat, gc, repeats]
    outputs = [output]
    options = {
        'memory': '10g',
        'walltime': '1-00:00:00',
        'account': 'CancerEvolution'
    }
    
    spec = f'''
    
    python create_feature_file.py {ascat} {centromere} {gc} {repeats} {output}
    
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def discretize_featurefile(featurefile, output, first_output): 
    inputs = [featurefile]
    outputs = [first_output]
    options = {
        'memory': '10g',
        'walltime': '1-00:00:00',
        'account': 'CancerEvolution'
    }
    
    spec = f'''
    
    python discretize.py {featurefile} {output}
    
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def gensimLDA(features, ntopics, bins):
    output = f'../steps/gensim/lda/lda_t{ntopics}_f{nfeat}_b{bins}.model'

    inputs = [features]
    outputs = [output]
    options = {
        'memory': '10g',
        'walltime': '1-00:00:00',
        'account': 'CancerEvolution'
    }

    spec = f'''
    
    python gensimLDA.py {features} {ntopics} {output}

    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def gensimHDP(features, nfeat):
    outputname = f'../steps/gensim/hdp/hdp_f{nfeat}.model'

    inputs = [features]
    outputs = [outputname]
    options = {
        'memory': '10g',
        'walltime': '1-00:00:00',
        'account': 'CancerEvolution'
    }

    spec = f'''
    
    python gensimHDP.py {features} {outputname}

    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def gensimNMF(features, ntopics, nfeat):
    outputname = f'../steps/gensim//nmf/nmf_t{ntopics}_f{nfeat}.model'

    inputs = [features]
    outputs = [outputname]
    options = {
        'memory': '10g',
        'walltime': '3-00:00:00',
        'account': 'CancerEvolution'
    }

    spec = f'''
    
    python gensimNMF.py {features} {ntopics} {outputname}

    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

gwf.target_from_template(
    name='UpdateAscat',
    template=update_ascat(
        samplefiles=sample,
        ascat=ascat,
        updatedascat=updatedascat
    )
)

gwf.target_from_template(
    name='SampleAscat',
    template=sample_ascat(
        ascat=updatedascat,
        nsamples=nsamples,
        output=sampledascat
    )
)

gwf.target_from_template(
    name='SampleAscat30',
    template=sample_ascat(
        ascat=updatedascat,
        nsamples=0.3,
        output=remaining30ascat
    )
)


gwf.target_from_template(
    name=f'CreateFeatures',
    template=create_feature_file(
        ascat = sampledascat,
        centromere = centromereinfo,
        gc = gc,
        repeats = repeats,
        output = featurefile
    )
)

gwf.target_from_template(
    name=f'CreateFeaturesValidation',
    template=create_feature_file(
        ascat = remaining30ascat,
        centromere = centromereinfo,
        gc = gc,
        repeats = repeats,
        output = validation_featurefile
    )
)

gwf.target_from_template(
    name=f'DiscretizeFeatures',
    template=discretize_featurefile(
        featurefile = featurefile,
        output = discrete_features,
        first_output = f'{discrete_features}_8.features'
    )
)

gwf.target_from_template(
    name=f'DiscretizeFeatures30',
    template=discretize_featurefile(
        featurefile = validation_featurefile,
        output = '../steps/featurefiles/discretized_10_0.3',
        first_output = '../steps/featurefiles/discretized_10_0.3_8.features'
    )
)

for i in range(start, ntopics + 1):

    gwf.target_from_template(
        name=f'gensimLDA_t{i}_f{nfeat}',
        template=gensimLDA(
            features=f'{discrete_features}_8.features',
            ntopics=i,
            bins=8
        )
    )

#     gwf.target_from_template(
#         name=f'gensimNMF_t{i}_f{nfeat}_{i}',
#         template=gensimNMF(
#             features=discrete_features,
#             ntopics=i,
#             nfeat=nfeat
#         )
#     )

# gwf.target_from_template(
#     name=f'gensimHDP_f{nfeat}_{i}',
#     template=gensimHDP(
#         features=discrete_features,
#         nfeat=nfeat
#     )
# )