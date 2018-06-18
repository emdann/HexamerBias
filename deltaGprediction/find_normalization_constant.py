import pandas as pd
import numpy as np
from hexVSprimed import *

ptMatrixCele="CG-pbat-gDNA-CeleTotal-noBS-1preAmp-handMix_lmerged_R1.ptCounts.qualFilt.csv.gz"
ptMatCele=pd.read_csv(ptMatrixCele, compression=findCompr(ptMatrix), index_col=0)
ptMatrixZfish="CG-pbat-gDNA-zfishTotal-noBS-1preAmp-handMix_lmerged_R1.ptCounts.qualFilt.csv.gz"
ptMatZfish=pd.read_csv(ptMatrixZfish, compression=findCompr(ptMatrix), index_col=0)
