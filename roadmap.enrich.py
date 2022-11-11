# Written by Yuwen Pan

import argparse
import pandas
import numpy as np 
import scipy.stats as stats
from concurrent.futures import ProcessPoolExecutor

## input
parser = argparse.ArgumentParser()
parser.add_argument("--cis", type=str, \
                    help="SNPs cis-regulating gene expression. 2 columns, <chromosome ID [e.g. chr1]> <position>, other columns will be ignored, no header.", \
                    required=True)
parser.add_argument("--bg", type=str, \
                    help="All SNPs in cis-region. The same input format.", \
                    required=True)
parser.add_argument("--out", type=str, \
                    help="/path/of/output/file/prefix", \
                    required=True)
args = parser.parse_args()

cisfile = args.cis
bgfile = args.bg
output = args.out

## annotate
globals()['annotationsource'] = pandas.read_csv('/picb/humpopg-bigdata5/panyuwen/Annotation/roadmap/15states5markers/blood/merged.annotation.bed.gz',sep='\t',header=None)

globals()['bgdata'] = pandas.read_csv(bgfile,sep='\s+',header=None)
globals()['bgdata'].index = globals()['bgdata'][0] +":"+ globals()['bgdata'][1].apply(str)
globals()['cisdata'] = pandas.read_csv(cisfile,sep='\s+',header=None)
globals()['cisdata'].index = globals()['cisdata'][0] +":"+ globals()['cisdata'][1].apply(str)

def bg_annotate(chromid):
    annotation = globals()['annotationsource'][globals()['annotationsource'][0]==chromid].copy()
    bgtemp = globals()['bgdata'][globals()['bgdata'][0]==chromid].copy()    
    bgtemp['anno'] = bgtemp.apply(lambda x: ';'.join(list(annotation[(annotation[1]<x[1]) & (annotation[2]>=x[1])][3].unique())),axis=1)
    return bgtemp

with ProcessPoolExecutor(max_workers = 22) as pool:
    #bg = list(pool.map(bg_annotate,['chr'+str(x) for x in range(1,23)]))
    bg = pandas.concat(list(pool.map(bg_annotate,list(set(bgdata[0])))))

cis = bg.loc[list(cisdata.index)]

bg[[0,1,'anno']].to_csv(output+'.background.roadmap.anno.txt',header=None,index=None,sep='\t')
cis[[0,1,'anno']].to_csv(output+'.cis.roadmap.anno.txt',header=None,index=None,sep='\t')

## count
states = ['1_TssA','2_TssAFlnk','3_TxFlnk','4_Tx','5_TxWk','6_EnhG','7_Enh','8_ZNF/Rpts','9_Het','10_TssBiv','11_BivFlnk','12_EnhBiv','13_ReprPC','14_ReprPCWk','15_Quies']
count = pandas.DataFrame(index=states, columns=['background','significant'])

for i in states:
    count.loc[i,'background'] = bg[bg['anno'].apply(lambda x: i in x)].shape[0]
    count.loc[i,'significant'] = cis[cis['anno'].apply(lambda x: i in x)].shape[0]
count.loc['total_sites',:] = [bg.shape[0],cis.shape[0]]

with open(output+'.annotation_sites_count.txt','w') as f:
    f.write('background: '+str(bg.shape[0])+' sites\n')
    f.write('significant: '+str(cis.shape[0])+' sites\n')

######################################################################
#            state   nonstate
#     eQTL     A         B     A+B     A = count.loc[state,'significant'];
#                                      B = count.loc['total_sites','significant'] - count.loc[state,'significant'];
# non-eQTL     C         D     C+D     C = count.loc[state,'background'] - count.loc[state,'significant'];
#                                      D = count.loc['total_sites','background'] - count.loc['total_sites','significant'] - count.loc[state,'background'] + count.loc[state,'significant']
#             A+C       B+D
######################################################################
## fisher_exact([A,B],[C,D])
for state in states:
    A = count.loc[state,'significant']
    B = count.loc['total_sites','significant'] - count.loc[state,'significant']
    C = count.loc[state,'background'] - count.loc[state,'significant']
    D = count.loc['total_sites','background'] - count.loc['total_sites','significant'] - count.loc[state,'background'] + count.loc[state,'significant']
    count.loc[state,'OR'], count.loc[state,'pvalue'] = stats.fisher_exact([[A,B],[C,D]])
count.drop('total_sites',inplace=True)

count.to_csv(output+'.RoadMap_enrich_stat.txt',sep='\t')

## for plot
statesfullname = ['Active TSS','Flanking Active TSS',"Transcr. at gene 5' and 3'",'Strong transcription','Weak transcription','Genic enhancers','Enhancers','ZNF genes & repeats','Heterochromatin','Bivalent/Poised TSS','Flanking Bivalent TSS/Enh','Bivalent Enhancer','Repressed PolyComb','Weak Repressed PolyComb','Quiescent/Low']
count['state'] = statesfullname
count = count[['OR','pvalue','state']]
count['log10pvalue'] = -1.0 * np.log10(np.array(count['pvalue']))
count.to_csv(output+'.RoadMap_enrich_stat.4plot.txt',sep='\t',index=None)


