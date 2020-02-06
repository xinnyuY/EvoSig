import gzip,pickle
import os
from pandas import DataFrame
import matplotlib.pyplot as plt
import csv

###### Modification manually*
output = '/home/xinyu/Downloads/ED_project/TCGA/4.rank_result/2020-01-06/'

######NMF by type 
for filename in os.listdir(output):
  
  summary = DataFrame(columns=['rank', 'cophenetic', 'dispersion', 'evar', 'rss', 'sparseness1', 'sparseness2','euclidean', 'kl'])
  
  for f in os.listdir(output+filename+'/'):
   	with open(output+filename+'/'+f,'rb') as in_f:
	       rank = pickle.load(in_f)
	       index_rank = int(f.split('_')[2].split('.')[0])
	       rank_estimate = [rank[index_rank]['rank'],rank[index_rank]['cophenetic'],rank[index_rank]['dispersion'],rank[index_rank]['evar'],rank[index_rank]['rss'],rank[index_rank]['sparseness'][0],rank[index_rank]['sparseness'][1],rank[index_rank]['euclidean'],rank[index_rank]['kl']]
	       summary.loc[-1] =rank_estimate
	       summary.index = summary.index + 1
	       summary = summary.sort_index()
  # control output path 
  csv.register_dialect('rank_dia',skipinitialspace=True)
  ###### Modification manually*
  export_csv = summary.to_csv(r'/home/xinyu/Downloads/ED_project/TCGA/4.rank_result/2020-01-06/'+filename+'.csv', index = None, header=True)
