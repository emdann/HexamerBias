import os
import pandas as pd

dir='predictedDg/'
for file in os.listdir(dir):
    mat=pd.read_csv(dir+file, index_col=0)
    print(file.split('_')[1], "Max:",str(mat.values.max()), "Min:",str(mat.values.min()))
