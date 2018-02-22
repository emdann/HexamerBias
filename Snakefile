
SAMPLE = 'SvdB11d1-MitoTrackerThird-Satellites-Adult'
TYPE = 'rna'
DIR = '/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/testing/'
REFGEN='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked_eGFP_Mito.fa'
CELLS=['cell100', 'cell2']
# rule all:
#     input:
#         # html=expand("fastqc/{sample}_{read}_fastqc.html",sample=SAMPLE, read=READS),
#         met_extraction=expand("met_extraction/{sample}_{read}_bismark_bt2.deduplicated.CpG_report.txt.gz", sample=SAMPLE, read=READS)
#         # split_bam=expand("bam/{sample}_{read}_bismark_bt2.bam", sample=SAMPLE, read=READS)
#         # bam1=expand("bam/{sample}_1_bismark_bt2.bam", sample=SAMPLE)


rule get_primed_region:
    input:
        bam='/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/{sample}.sam.gz'
    output:
        primedfa= '{dir}/{sample}.primedreg.fa'
    params:
        refgen=REFGEN,
        t=TYPE,
        out=DIR
    script:
        "getPrimedRegion.py -o {params.out} {input.bam} {params.refgen} {params.t}"

rule kmer_count:
    input:
        coutt='/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/{sample}.coutt.csv',
        refgen=REFGEN
    output:
        cellAbundance='{dir}/{sample}.cellAbundance.noN.csv'
    params:
        out={dir}
    threads: 8
    script:
        "kmerCounter.py -o {params.out} {input.coutt} {input.refgen}"

rule pt_counts:
    input:
        bam='/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/{sample}.sam.gz',
        primedfa= '{dir}/{sample}.primedreg.fa'
    output:
        ptCounts=expand('{dir}/ptCounts/{sample}.{cell}.ptCounts.qualFilt.parallel.csv', cell=CELLS)
    threads: 10
    script:
        "cellPrimerTemplTab.py -o {input.bam} {input.primedfa}"
