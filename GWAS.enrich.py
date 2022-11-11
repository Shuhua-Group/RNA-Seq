# Written by Ke Huang

import pandas as pd
import numpy as np
import os,sys

'''transform pvalue to float'''

def GWAS_p(cutoff,out):
    pvalue_list=list()
    data=pd.read_table('/picb/humpopg-bigdata5/huangke/research/Sex_bias/XJU90_Sex_eQTL/55.EAS.EUR.biased.adm.mediate/real/GWAS/GWAS_info.txt')
    number1=data['P-VALUE'].str.split('-').str[0].str[0]
    number2=data['P-VALUE'].str.split('-').str[1]
    for i in range(len(data)):
        try:
            pvalue=int(number1[i])*pow(10,-int(number2[i]))
            pvalue_list.append(pvalue) 
        except:
            print(data.iloc[i])
            pvalue_list.append(0)
    data['pvalue']=pvalue_list
    #data=data.loc[data['pvalue']<pow(10,-8)]
    data=data.loc[data['pvalue']<cutoff]
    print(len(data['SNP_GENE_IDS'].drop_duplicates()))
    
    ## overlap gene
    data=data.reset_index(drop=True)
    df_list=list()
    for j in range(len(data)):
        try:
            gene_list=data['SNP_GENE_IDS'].iloc[j].split(',')
            if len(gene_list)>1:
                df=pd.DataFrame(data.iloc[j]).T
                
                data=data.drop([j])
                for gene in gene_list:
                    df['SNP_GENE_IDS']=gene
                    df_list.append(df)
        except:
            pass
    data=pd.concat([pd.concat(df_list),data])

    data.to_csv('./cutoff'+str(out)+'/GWAS.info.pvalue.txt',index=None,sep='\t')

''' GWAS gene overlap with sexpop gene'''

def GWAS_gene(out,filename):
    data_merge=pd.read_table(filename)
    data_merge['SNP_GENE_IDS']=data_merge['geneid']
    print(len(data_merge))
    data_GWAS=pd.read_table('./cutoff'+str(out)+'/GWAS.info.pvalue.txt')
    data_merge=pd.merge(data_GWAS,data_merge)
    data_merge.to_csv('./cutoff'+str(out)+'/GWAS.accorand.model.gene.txt',index=None,sep='\t')

''' enrich number table'''
def enrich_number(out):
    data_GWAS=pd.read_table('./cutoff'+str(out)+'/GWAS.info.pvalue.txt')
    df=pd.DataFrame(data_GWAS['DISEASE/TRAIT'].value_counts())
    df['GWAS_number']=df['DISEASE/TRAIT']
    df['DISEASE/TRAIT']=df.index
    df['GWAS_total']=len(data_GWAS['SNP_GENE_IDS'].drop_duplicates())
    df['pvalue_cutoff']='10-8'
    df.to_csv('./cutoff'+str(out)+'/GWAS.info.dataset.txt',index=None,sep='\t')

    data_sexpop=pd.read_table('./cutoff'+str(out)+'/GWAS.accorand.model.gene.txt')
    df=pd.DataFrame(data_sexpop['DISEASE/TRAIT'].value_counts())
    df['loci_number']=df['DISEASE/TRAIT']
    df['DISEASE/TRAIT']=df.index
    df['loci_total']=len(data_sexpop['SNP_GENE_IDS'].drop_duplicates())
    df.to_csv('./cutoff'+str(out)+'/sexpop.info.dataset.txt',index=None,sep='\t')

''' enrich socre'''
def fold_enrich(out):
    import scipy.stats as stats
    odd_ratio_list=list()
    fisher_pvalue_list=list()
    data_GWAS=pd.read_table('./cutoff'+str(out)+'/GWAS.info.dataset.txt')
    data_loci=pd.read_table('./cutoff'+str(out)+'/sexpop.info.dataset.txt')
    data_merge=pd.merge(data_GWAS,data_loci)
    data_merge=data_merge.reset_index(drop=True)
    for i in range(len(data_merge)):
        number_list=np.array([[data_merge['loci_number'][i],data_merge['loci_total'][i]],[data_merge['GWAS_number'][i],data_merge['GWAS_total'][i]]])
        odd,pvalue=stats.fisher_exact(number_list)
        odd_ratio_list.append(odd)
        fisher_pvalue_list.append(pvalue)
    data_merge['odds']=odd_ratio_list
    data_merge['fisher.pvalue']=fisher_pvalue_list
    data_merge.to_csv('./cutoff'+str(out)+'/GWAS.loci.enrich.txt',index=None,sep='\t')
    data_merge=data_merge.loc[(data_merge['odds']>1)&(data_merge['fisher.pvalue']<0.05)&(data_merge['GWAS_number']>10)&(data_merge['GWAS_number']<500)]
    data_merge['-log10(GWAS cutoff pvalue)']=out

    data_merge.to_csv('./cutoff'+str(out)+'/GWAS.loci.enrich.sig.txt',index=None,sep='\t')


filename = sys.argv[1]
for cut in [8,10,15,20,25,30]:
    os.system('mkdir -p cutoff'+str(cut))
    GWAS_p(pow(10,-cut),cut)
    GWAS_gene(cut,filename)
    enrich_number(cut)

    fold_enrich(cut)

