
SAMPLE = 'SvdB11d1-MitoTrackerThird-Satellites-Adult'
TYPE = 'rna'
DIR = '/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/testing/'
REFGEN='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked_eGFP_Mito.fa'

rule get_primed_region:
    input:
        bam=expand('{{sample}}.bam')
    output:
        txt='{{sample}}.primedreg.fa'
    params:
        refgen=REF_GEN,
        t=TYPE,
        out=DIR
    script:
        "getPrimedRegion.py -o {params.out} {input.bam} {params.refgen} {params.t}"
