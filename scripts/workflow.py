from gwf import Workflow, AnonymousTarget
import modpath

gwf = Workflow()

"""
Workflow for LDA analysis of ASCAT data
"""

sample = '/home/janneae/TCGA/DerivedData/PanCancer/TCGA_ASCAT_RAW_PVL/ASCAT_TCGA/TCGA*.segments.raw.txt' 
centromereinfo = '../data/chrominfo.snp6.txt'
ascat = '../data/filteredAscat.txt'
gc = '../data/gc_content.txt'
output = '../steps/AllTCGA_features_mateo10.txt'

def create_feature_file(sample, centromere, ascat, gc, output): 
    inputs = [centromere, ascat, gc]
    outputs = [output]
    options = {
        'memory': '5g',
        'walltime': '05:00:00'
    }
    
    spec = f'''

    python create_feature_file.py {sample} {centromere} {ascat} {gc} {output}
    
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

gwf.target_from_template(
    name='create_features',
    template=create_feature_file(
        sample = sample,
        centromere = centromereinfo,
        ascat = ascat,
        gc = gc,
        output = output
    )
)