from subsmpCsv import rand_subsample_csv as subrand
import pandas as pd

zf_file='/hpc/hub_oudenaarden/edann/hexamers/zf_prediction/zebrafishBs.commonPtPairs.csv'
mm_file = '/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/VAN1667.commonPtPairs.csv'

# df_zf = rand_subsample_csv(zf_file, index_col=0)
# df_mm = pd.read_csv(mm_file, index_col=0)
# subsmp_df_mm = df_mm.loc[[i for i in df_zf.index if i in df_mm.index]]
# c = pd.concat([subsmp_df_mm, df_zf], axis =1)

## Build correlation matrix
corrMatrix = {}
n = sum(1 for line in open(filename)) - 1 # get number of lines
for ncol_zf in range(len(open(zf_file).readline().strip().split(',')))[1:]:
    df_zf = pd.read_csv(zf_file, usecols=[0,ncol_zf], index_col=0)
    zf_sample = df_zf.columns[0]
    print("Loaded sample %s" % zf_sample)
    print()
    corrMatrix[zf_sample]={}
    for ncol_mm in range(len(open(mm_file).readline().strip().split(',')))[1:]:
        df_mm = pd.read_csv(mm_file, usecols=[0,ncol_mm], index_col=0)
        mm_sample = df_mm.columns[0]
        c = pd.concat([df_mm, df_zf], axis =1)
        pcc = c.corr().iloc[0,1]
        corrMatrix[zf_sample][mm_sample] = pcc

pd.DataFrame(corrMatrix).to_csv('correlationMatrix.zfMm.csv')
