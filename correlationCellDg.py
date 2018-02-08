import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from scipy.stats.stats import pearsonr as pcc

mat1file='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/cell121_ptDg.csv'
mat2file='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/cell147_ptDg.csv'

mat1file='~/mnt/edann/hexamers/rnaseq/cell121_ptDg.csv'
mat2file='~/mnt/edann/hexamers/rnaseq/cell147_ptDg.csv'


mat1 = pd.read_csv(mat1file, sep=',', index_col=0)
mat2 = pd.read_csv(mat2file, sep=',', index_col=0)

array1 = np.array(mat1)
array2 = np.array(mat2)

xAx=[]
yAx=[]
for (x,y), value in np.ndenumerate(array1):
    # print(value, array2[x][y])
    if value != -99999.0 and array2[x][y] != -99999.0:
        xAx.append(value)
        yAx.append(array2[x][y])

plt.figure(4)
plt.plot(xAx, yAx, 'b.')
plt.xlabel('deltaG cell 121')
plt.ylabel('deltaG cell 147')
plt.savefig('/hpc/hub_oudenaarden/edann/output/test_pyplot_cell121_cell147.pdf')
