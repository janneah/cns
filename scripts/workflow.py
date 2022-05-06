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
featurefile = f'../steps/featurefiles/discretized_{nfeat}_{nsamples}_4.features'

# Validation input and output
remaining30ascat = '../steps/sampledascat/sampled_0.3.ascat'
validation_featurefile = f'../steps/featurefiles/discretized_{nfeat}_0.3_6.features'

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

def gensimLDA(features, ntopics, nfeat):
    outputname = f'../steps/gensim/lda/lda_t{ntopics}_f{nfeat}_8bins.model'

    inputs = [features]
    outputs = [outputname]
    options = {
        'memory': '10g',
        'walltime': '1-00:00:00',
        'account': 'CancerEvolution'
    }

    spec = f'''
    
    python gensimLDA.py {features} {ntopics} {outputname}

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
    name=f'CreateFeatures_4',
    template=create_feature_file(
        ascat = sampledascat,
        centromere = centromereinfo,
        gc = gc,
        repeats = repeats,
        output = featurefile
    )
)

gwf.target_from_template(
    name='CreateFeaturesValidation',
    template=create_feature_file(
        ascat = remaining30ascat,
        centromere = centromereinfo,
        gc = gc,
        repeats = repeats,
        output = validation_featurefile
    )
)

for i in range(start, ntopics + 1):

    gwf.target_from_template(
        name=f'gensimLDA_t{i}_f{nfeat}',
        template=gensimLDA(
            features=featurefile,
            ntopics=i,
            nfeat=nfeat
        )
    )

    gwf.target_from_template(
        name=f'gensimNMF_t{i}_f{nfeat}',
        template=gensimNMF(
            features=featurefile,
            ntopics=i,
            nfeat=nfeat
        )
    )

gwf.target_from_template(
    name=f'gensimHDP_f{nfeat}',
    template=gensimHDP(
        features=featurefile,
        nfeat=nfeat
    )
)