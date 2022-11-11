# Written by Yuwen Pan

import argparse
import pandas
import numpy as np
import scipy
import scipy.stats as stats
from math import factorial
from math import exp
from scipy.stats import poisson
# Written by Yuwen Pan

from concurrent.futures import ProcessPoolExecutor

## input
parser = argparse.ArgumentParser()
parser.add_argument("--cis", type=str, \
					help="SNPs cis-regulating gene expression. 2 columns, <chromosome ID [e.g. chr1]> <position>, other columns will be ignored, no header.", \
					required=True)
parser.add_argument('--window_overlap_size',type=str, \
					help="window size and overlap length, format: window.size@overlap.length, e.g. 1000000@100000.", \
					required=False, default='1000000@100000')
parser.add_argument("--out", type=str, \
					help="/path/of/output/file/prefix", \
					required=True)
args = parser.parse_args()

cisfile = args.cis
window_overlap = args.window_overlap_size
output = args.out

## main
size = pandas.read_csv('/picb/humpopg-bigdata5/tanxinjiang/QTL_Annotation/chrom.size.hg19.ucsc.txt',sep='\t',header=None,names=['chr','len'],index_col='chr')
windowsize, overlapsize = window_overlap.strip().split('@')
globals()['windowsize'] = int(windowsize); globals()['overlapsize'] = int(overlapsize); globals()['stepsize'] = windowsize-overlapsize

def counteqtl(chromid):
	temp = pandas.DataFrame(columns=['chr','start','end','len','eqtl','pvalue','adj.pvalue'])
	cis = pandas.read_csv(cisfile,sep='\s+',header=None,names=['chr','pos'])
	loc = cis[cis['chr'] == 'chr'+str(chromid)].copy()
	if loc.shape[0] == 0:
		return temp
	length = globals()['windowsize']
	step = globals()['stepsize']
	temp['start'] = range(0,size.loc['chr'+str(chromid),'len'],step)
	temp['end'] = temp['start'] + length
	if temp.iloc[temp.shape[0]-1,2] > size.loc['chr'+str(chromid),'len']:
		temp.iloc[temp.shape[0]-1,2] = size.loc['chr'+str(chromid),'len']
	temp['len'] = temp['end'] - temp['start']
	temp['chr'] = 'chr' + str(chromid)
	temp['eqtl'] = temp.apply(lambda x: loc[(loc['pos']>=x['start']) & (loc['pos']<x['end'])].shape[0],axis=1)
	## Poisson distribution: exp(-1.0 * lambda) * pow(lambda,K) / factorial(K)
	total = loc.shape[0]
	for i in range(temp.shape[0]):
		lam = temp.loc[i,'len'] * 1.0/size.loc['chr'+str(chromid),'len'] * total
		temp.loc[i,'pvalue'] = scipy.stats.distributions.poisson.pmf(temp.loc[i,'eqtl'], lam)
	return temp
	
with ProcessPoolExecutor(max_workers = 23) as pool:
	count = pandas.concat(list(pool.map(counteqtl,range(1,24))))
count.index = count['chr'].apply(str) +":"+ count['start'].apply(str)

def fdr(pvalues):
	## the pvalue should have already been sorted
	## and pvalues should be pandas.Series
	n = len(pvalues)
	pvalues.dropna(inplace=True)
	adj_pvalues = pvalues * len(pvalues) / range(1,len(pvalues)+1)
	if adj_pvalues[-1] > 1.0:
		adj_pvalues[-1] = 1.0
	for i in range(len(adj_pvalues)-2,-1,-1):
		adj_pvalues[i] = min(adj_pvalues[i+1],adj_pvalues[i])
	adj_pvalues = list(adj_pvalues) + [np.nan] * (n-len(adj_pvalues))
	return adj_pvalues

count.sort_values(by='pvalue',ascending=True,inplace=True)
count['adj.pvalue'] = np.array(fdr(count['pvalue']))
count.sort_values(by=['chr','start'],ascending=True,inplace=True)
count['chr'] = count['chr'].apply(lambda x: x[3:])

count.to_csv(output,sep='\t',index=None)
