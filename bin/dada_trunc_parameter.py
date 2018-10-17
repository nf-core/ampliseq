#!/usr/bin/env python3
#@author Alexander Peltzer
# Takes multiple CSV files and generates a tuple of index locations > that resemble the cutoff value used for DADA2 in QIIME2.

import pandas as pd
import sys 


data_fw = pd.read_csv(sys.argv[1]) 
data_rv = pd.read_csv(sys.argv[2])

