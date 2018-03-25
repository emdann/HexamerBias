import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

def from_refgen_to_bed(refgen):
    '''
    Make bed like entries for chromosomes in reference genome fasta file
    '''
    beds=[]
    with ps.FastxFile(refgen) as f:
        for entry in f:
            beds.append((entry.name + ' 0 ' + str(len(entry.sequence))))
    return(beds)

def genome_wide_artificial_coverage(covFile,refgen,abundanceFile, outfile):
    abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
    coverage = pd.read_table(covFile, index_col=0, sep='\t',compression=findCompr(covFile))
    density = template_density(coverage.exp,abundance)
    beds = from_refgen_to_bed(refgen)
    save_bigWig(beds,refgen_fasta, outfile = outfile, threads=10)
    return('')

abundanceFile = "/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz"
covFile = "predictedCoverage_avgVAN1667.txt"
refgen='/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa'
genome_wide_artificial_coverage(covFile, refgen, abundanceFile, outfile='mm10.artCov.bw')
#
# bedFile="/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.random.40.bed"
