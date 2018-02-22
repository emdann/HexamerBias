
SAMPLE = 'SvdB11d1-MitoTrackerThird-Satellites-Adult'
TYPE = 'rna'
DIR = '/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/testing/'
REFGEN='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked_eGFP_Mito.fa'

rule get_primed_region:
    input:
        bam='/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/{sample}.sam.gz'
    output:
        primedfa=expand(DIR+'{sample}.primedreg.fa', dir=DIR)
    params:
        refgen=REFGEN,
        t=TYPE,
        out=DIR
    script:
        "getPrimedRegion.py -o {params.out} {input.bam} {params.refgen} {params.t}"
