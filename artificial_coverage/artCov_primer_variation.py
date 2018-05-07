import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/optimization')
from primerProbability import *
from predictCovBs import *
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

def progressive_variation(ppm, step=0.05, increase = True, nuc='A',pos=0):
    '''
    Makes ppm with progressively higher concentration of one base at one position
    '''
    ppm1 = ppm.copy()
    if increase:
        ppm1[pos][nuc]+=step
    else:
        ppm1[pos][nuc]-=step
    ppm1[pos][nuc] = ppm1[pos][nuc].round(2)
    newProbs = (1-ppm1[pos][nuc])/3
    ppm1[pos][ppm1[pos].index!=nuc] = newProbs
    return(ppm1)

def build_variation_probabilities(ppm, step=0.05, increase = True, nuc='A',pos=0):
    varDf = get_even_prob()
    if increase:
        while ppm[pos][nuc] < 1:
            ppm = progressive_variation(ppm, step=step, increase = increase, nuc=nuc,pos=pos)
            newProb = prob_from_ppm(ppm, all_hexamers())
            newProb.columns = [nuc + '_frac' + str(ppm[pos][nuc])]
            varDf = pd.concat([varDf,newProb], axis=1)
    else:
        while ppm[pos][nuc] > 0:
            ppm = progressive_variation(ppm, step=step, increase = increase, nuc=nuc,pos=pos)
            newProb = prob_from_ppm(ppm, all_hexamers())
            newProb.columns = [nuc + '_frac' + str(ppm[pos][nuc])]
            varDf = pd.concat([varDf,newProb], axis=1)
    return(varDf)

# def artificial_cov_bed_entry(params):
#     bedEntry,refgen_fasta,density,read_length=params
#     chr,start,end = bedEntry.split()
#     print('Processing entry ', bedEntry, flush=True)
#     seq = ps.FastaFile(refgen_fasta).fetch(reference=chr, start=int(start), end=int(end)).upper()
#     strandSpecificPosDic = sum_strands_per_base_cov(seq, density, int(start), readLength=read_length)
#     smoothPosDic = kernel_smoothing(strandSpecificPosDic)
#     return(chr,smoothPosDic)

def make_multi_covprofile_from_cov_df(varDf, bedEntry, refgen, abundance):
    artCovDf = pd.DataFrame()
    for cov in varDf.iteritems():
        name,coverage=cov
        density = template_density(coverage,abundance)
        artCov = artificial_cov_bed_entry((bedEntry,refgen,density,75))
        df = pd.DataFrame(artCov[1], index=[1]).T
        df.columns = [name]
        artCovDf = pd.concat([artCovDf, df], axis=1)
    return(artCovDf)

x = np.array([[0.25,0.25,0.25, 0.25, 0.25, 0.25],
            [0.25,0.25,0.25, 0.25, 0.25, 0.25],
            [0.25,0.25,0.25, 0.25, 0.25, 0.25],
            [0.25,0.25,0.25, 0.25, 0.25, 0.25]  ])

ppm = from_vec_to_ppm(x)

abundanceFile="/hpc/hub_oudenaarden/edann/hexamers/genomes_kmers/mm10.kmerAbundance.csv"
bedFile='artificial_coverage/highcov.random.42.bed'
refgen='/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa'

abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
with open(bedFile, 'r') as f:
    beds = [line.strip() for line in f.readlines()]
bedEntry=beds[2]

varPos0 = pd.DataFrame()
for nuc in ['A', 'T', 'C', 'G']:
    varDf = build_variation_probabilities(ppm, nuc=nuc, step=0.1)
    nucDf = make_multi_covprofile_from_cov_df(varDf, bedEntry, refgen, abundance)
    varPos0 = pd.concat([varPos0, nucDf], axis=1)
