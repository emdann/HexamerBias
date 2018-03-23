import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

def from_refgen_to_bed(genomefile):
    '''
    Make bed like entries for chromosomes in genome file
    '''
    beds=[]
    with open(genomefile, 'r') as f:
        for line in f.readlines():
            chr,length = line.strip().split()
            if '_' not in chr:
                beds.append(chr+' 0 '+length)
    return(beds)

def genome_wide_artificial_coverage(genomefile,covFile,refgen,abundanceFile, outfile):
abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
coverage = pd.read_table(covFile, index_col=0, sep='\t',compression=findCompr(covFile))
density = template_density(coverage.exp,abundance)
beds = from_refgen_to_bed(genomefile)
covBed = artificial_cov(beds,refgen,density)
save_coverage_bed(covBed, outfile=outfile)
    return(covBed)

genomefile='/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.genome'
abundanceFile = "/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz"
covFile = "predictedCoverage_avgVAN1667.txt"
refgen='/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa'
genome_wide_artificial_coverage(genomefile, covFile, refgen,abundanceFile, outfile='mm10.artificial_cov.VAN1667avg.bed')

bedFile="/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.random.40.bed"
