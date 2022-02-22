from gwf import Workflow, AnonymousTarget

gwf = Workflow()

"""
Workflow for LDA analysis of ASCAT data
"""

sample = '/home/janneae/TCGA/DerivedData/PanCancer/TCGA_ASCAT_RAW_PVL/ASCAT_TCGA/TCGA*.segments.raw.txt' 
centromereinfo = '../data/chrominfo.snp6.txt'
ascat = '../data/filteredAscat.txt'
gc = '../data/gc.content.txt'
output = '../steps/AllTCGA_features_mateo10.txt'
updatedascat = '../data/filteredAscatRaw.txt'


def update_ascat(samplefiles, ascat, output):
    inputs = [ascat]
    outputs = [output]
    options = {
        'memory': '8g',
        'walltime': '24:00:00'
    }
    spec = f'''
    
    python update_ascat.py {samplefiles} {ascat} {output}

    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def create_feature_file(ascat, centromere, gc, output): 
    inputs = [centromere, ascat, gc]
    outputs = [output]
    options = {
        'memory': '5g',
        'walltime': '05:00:00'
    }
    
    spec = f'''

    python create_feature_file.py {ascat} {centromere} {gc} {output}
    
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

gwf.target_from_template(
    name = 'UpdateAscat',
    template = update_ascat(
        samplefiles=sample,
        ascat=ascat,
        output=updatedascat
    )
)

gwf.target_from_template(
    name='CreateFeatures',
    template=create_feature_file(
        ascat = updatedascat,
        centromere = centromereinfo,
        gc = gc,
        output = output
    )
)