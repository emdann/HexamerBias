import pandas as pd
from hexVSprimed import *

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Count pt events in ptMatrix \n By Emma Dann")
argparser.add_argument('ptmatrix', type=str, help='ptCounts file')
args = argparser.parse_args()

ptMatrix=args.ptmatrix

ptMat = pd.read_csv(ptMatrix, compression=findCompr(ptMatrix), index_col=0)
print(ptMat.sum().sum())
