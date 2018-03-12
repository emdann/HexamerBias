#!/bin/sh
# properties = {"rule": "get_primed_region", "local": false, "input": ["/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/SvdB11d2-MitoTrackerThird-Satellites-Adult.sam.gz"], "output": ["/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/test_snakemake/SvdB11d2-MitoTrackerThird-Satellites-Adult.primedreg.fa"], "wildcards": ["/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/test_snakemake", "SvdB11d2-MitoTrackerThird-Satellites-Adult"], "params": {"refgen": "/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked_eGFP_Mito.fa", "t": "rna", "out": "/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/test_snakemake"}, "log": [], "threads": 1, "resources": {}, "jobid": 388, "cluster": {}}
cd /hpc/hub_oudenaarden/edann/bin/coverage_bias && \
/hpc/hub_oudenaarden/edann/venv2/bin/python3 \
-m snakemake /hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/test_snakemake/SvdB11d2-MitoTrackerThird-Satellites-Adult.primedreg.fa --snakefile /hpc/hub_oudenaarden/edann/bin/coverage_bias/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /hpc/hub_oudenaarden/edann/bin/coverage_bias/.snakemake/tmp.bviysqxs /hpc/hub_oudenaarden/aalemany/emma-adi/mouse/SvdB11d2-MitoTrackerThird-Satellites-Adult.sam.gz --latency-wait 5 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp  --force-use-threads  --allowed-rules get_primed_region  && touch "/hpc/hub_oudenaarden/edann/bin/coverage_bias/.snakemake/tmp.bviysqxs/388.jobfinished" || (touch "/hpc/hub_oudenaarden/edann/bin/coverage_bias/.snakemake/tmp.bviysqxs/388.jobfailed"; exit 1)

