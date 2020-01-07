import numpy as np
import pandas as pd
import nimfa
import os
from LoadOrWriteObj import write_or_read_obj


def read_rdata_csv(file_pa):

    raw_data = pd.read_csv(file_pa)
    raw_data = np.array(raw_data)
    raw_data = [list(row) for row in list(raw_data)]
    for row in raw_data:
        row.pop(0)
    data = np.array(raw_data)
    return data


#############test for per type ##############
path = '/home/xinyu/Downloads/ED_project/TCGA/3.ccfMatrix/2020-01-06/rank_estimate/'
output = '/home/xinyu/Downloads/ED_project/TCGA/4.rank_result/2020-01-06/'

for filename in os.listdir(path):
   data = read_rdata_csv(path+filename)
   cancertype =  filename.split("_")[0]
   formt = filename.split("_")[3]
   try:
     t_nmf = nimfa.Nmf(data)
     if os.path.isdir(output) == False : os.mkdir(output)
     if os.path.isdir(output+cancertype+"_"+formt) == False : os.mkdir(output+cancertype+"_"+formt) 
     else: continue
     for i in range(2, 19):
       rst_rank = t_nmf.estimate_rank(rank_range=[i], n_run=1000, what='all')
       write_or_read_obj(output+cancertype+"_"+formt+'/rst_for_'+str(i), rst_rank)
   except (RuntimeError, TypeError, NameError):
     pass
   
