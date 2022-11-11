#coding=utf-8

ind_info=open("ind_list_HU.txt")
ind_list=[]
for line in ind_info:
    line_s=line.split()
    ind_list.append(line_s[0])

file_out=open("expr_merge_unstrand.txt","w")
file_out.writelines("%s\t"%("gene_id"))
for i in range(len(ind_list)):
    file_out.writelines("%s\t"%ind_list[i])
file_out.writelines("\n")

gene_name=[]
file_tmp=open("/picb/humpopg-bigdata5/ningzhilin/RNA_Seq/RSEM_unstrand/WGC019716R/WGC019716R.genes.results")
file_tmp.readline()
for line in file_tmp:
    line_s=line.split()
    gene_name.append(line_s[0])

#exp_list=[gene1[],gene2[]]
exp_list=[]
for i in range(len(gene_name)):
    exp_list.append([])

for i in range(len(ind_list)):
    file="/picb/humpopg-bigdata5/ningzhilin/RNA_Seq/RSEM_unstrand/"+ind_list[i]+"/"+ind_list[i]+".genes.results" 
    f=open(file)
    f.readline()
    j=0
    for line in f:
        line_s=line.split()
        fpkm=line_s[-1]
        exp_list[j].append(fpkm)
        j+=1

for i in range(len(gene_name)):
    file_out.writelines("%s\t"%gene_name[i])
    for j in range(len(exp_list[i])):
        file_out.writelines("%s\t"%exp_list[i][j])
    file_out.writelines("\n")

file_out.close()

    
