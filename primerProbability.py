import numpy as np
import pandas as pd
'''
Compute the probability for all possible hexamers base on manufacturers concentrations.
Input: matrix of base concentration per position
Output: probability for all hexamers
'''

# Make structured array for concentrations
x = np.array([[0.25,0.25,0.25, 0.25, 0.25, 0.25],
            [0.25,0.25,0.25, 0.25, 0.25, 0.25]],
    dtype=[('foo', 'i4'),('bar', 'f4'), ('baz', 'S10')])
