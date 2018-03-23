#### ARTIFICIAL COVERAGE PROFILE
From the predicted coverage I want to make a igv like coverage track based on density of coverage for every hexamer (C/T).

Useful python functions are in module ```cov_from_density```.

## Making it work genome wide
Probably not a good idea to do the coverage of the whole genome. Nevertheless I am making this huge bed file with
```
source /hpc/hub_oudenaarden/edann/venv2/bin/activate; python /hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage/genome_wide_artificialcov.py
```



## Comparison of artificial and real coverage
