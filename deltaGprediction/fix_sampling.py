import pandas as pd
import numpy as np


def test_summing(conc, s):
    p = (conc.ab-conc.sumTemp+s)/sum(conc.ab)
    plog=p
    return(0.00001/max(plog), 0.00001/min(plog), max(plog)-min(plog))
