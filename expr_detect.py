# Written by Zhilin Ning

#coding=utf-8
from math import log

file_in=open("/picb/humpopg-bigdata5/ningzhilin/RNA_Seq/expr_merge.txt")
file_out=open("/picb/humpopg-bigdata5/ningzhilin/RNA_Seq/expr_detect.txt","w")
title=file_in.readline()
file_out.writelines("%s"%title)
for line in file_in:
    line_s=line.split()
    gene_name=line_s[0]
    for i in range(1,len(line_s)):
        line_s[i]=log((float(line_s[i])+1),2)
    mm=max(line_s[1:])
    if  mm > 1:
        for j in range(len(line_s)):
            file_out.writelines("%s\t"%line_s[j])
        file_out.writelines("\n")

file_out.close()

            
