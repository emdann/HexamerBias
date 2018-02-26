
SAMPLE = "SvdB11d2-MitoTrackerThird-Satellites-Adult"
TYPE = "rna"
DIR = "/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/test_snakemake"
REFGEN = "/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked_eGFP_Mito.fa"
CELLS = ['cell' + str(n) for n in range(1,385)]

rule all:
    input:
        predictedDg=expand('{dir}/predictedDg/{sample}_{cell}_ptDg_qual.csv', cell=CELLS, dir=DIR, sample=SAMPLE),
        numReads=expand('{dir}/{sample}.numReads.txt', sample=SAMPLE, dir=DIR)
        # split_bam=expand("bam/{sample}_{read}_bismark_bt2.bam", sample=SAMPLE, read=READS)
        # bam1=expand("bam/{sample}_1_bismark_bt2.bam", sample=SAMPLE)


rule get_primed_region:
    input:
        bam='/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/{sample}.sam.gz',
        refgen=REFGEN
    output:
        primedfa= '{dir}/{sample}.primedreg.fa'
    params:
        t=TYPE,
    #     out=DIR
    threads: 1
    script:
        "getPrimedRegion.py -o {wildcards.dir} {input.bam} {input.refgen} {params.t}"

rule kmer_count:
    input:
        coutt='/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/{sample}.coutt.csv',
        refgen=REFGEN
    output:
        cellAbundance='{dir}/{sample}.cellAbundance.noN.csv'
    # params:
    #     out=DIR
    threads: 8
    script:
        "kmerCounter.py -o {wildcards.dir} {input.coutt} {input.refgen}"

rule num_reads:
    input:
        bam='/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/{sample}.sam.gz'
    output:
        numReads='{dir}/{sample}.numReads.txt'
    # params:
    #     out=DIR
    threads: 1
    script:
        "numReadsPerCell.py -o {wildcards.dir} {input.bam}"

rule pt_counts:
    input:
        bam='/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/{sample}.sam.gz',
        primedfa= '{dir}/{sample}.primedreg.fa'
    output:
        ptCounts=expand('{{dir}}/ptCounts/{{sample}}.{cell}.ptCounts.qualFilt.parallel.csv', cell=CELLS)
    threads: 10
    script:
        "cellPrimerTemplTab.py -o {input.bam} {input.primedfa}"

rule predict_dg:
    input:
        ptCounts='{dir}/ptCounts/{sample}.{cell}.ptCounts.qualFilt.parallel.csv',
        cellAbundance='{dir}/{sample}.cellAbundance.noN.csv'
    output:
        predictedDg='{dir}/predictedDg/{sample}_{cell}_ptDg_qual.csv'
    threads: 1
    script:
        "ptModel.py {input.ptCounts} {input.cellAbundance}"
