## Prediction of Delta G from sequencing data
Scripts to compute predicted Delta G for hexamer binding energy from sequencing reads. For details of model see report.

First count kmers in reference genome with ```kmersInGenome.sh```, then ```run_deltaF_prediction``` runs the python scripts as described in workflow:

![workflow](DeltaGprediction_workflow.001.png)

Outputs matrix of predicted DeltaG for all possible template-primer pairs, assuming even concentration of hexamers.

**N.B. IN ALL MATRIXES TEMPLATES ARE IN ROWS AND PRIMERS ARE IN COLUMNS**  
