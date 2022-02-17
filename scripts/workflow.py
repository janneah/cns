from gwf import Workflow, AnonymousTarget
import glob

gwf = Workflow()

# Input files
sample = glob.glob('/home/janneae/cns/data/TCGA-*.segments.raw.txt') 
centromereinfo = '/home/janneae/cns/data/chrominfo.snp6.txt'
ascat = '/home/janneae/cns/data/filteredAscat.txt'
gc = '/home/janneae/cns/data/gc.content.txt'
output = '/home/janneae/cns/steps/TCGA_features_mateo_10.txt'

# '/home/janneae/TCGA/DerivedData/PanCancer/TCGA_ASCAT_RAW_PVL/ASCAT_TCGA/TCGA-*.segments.raw.txt'

def create_feature_file(sample, centromere, ascat, gc, output):    
    inputs = [sample, centromere, ascat, gc]
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