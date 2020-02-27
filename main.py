import csv
import numpy as np
import os
from bam_tools import *
from report_tools import *
bamfile = open("SRR8283737_cutviewbam.txt")
referencefile = open("Human IGHV F+ORF+in-frame P.fa")

reference_genelist = referencefile.readlines()
line = bamfile.readline()
line_list = line.split('\t')
if not os.path.exists('nucleotideCounts/'):
    os.makedirs('nucleotideCounts/')  # nucleotideCounts文件夹存放碱基信息
if not os.path.exists('MutCounts/'):
    os.makedirs('MutCounts/')  # MutCounts文件夹存放突变碱基信息


for i in range(0, len(reference_genelist)//2):  #2i是基因名，2i+1是基因
    fastqseq_list = []
    gene_header = reference_genelist[2*i].strip()
    reference_gene = reference_genelist[2*i+1].strip()
    nt_list = np.zeros((len(reference_gene), 4))  # nt_list用于存放碱基统计信息，参考基因组的长度行，4列，行索引代表碱基位置，列索引代表ATCG
    while True:
        if not line:
            break
        if line_list[2] not in gene_header:
            break

        fastqseq, quality = align(reference_gene, line_list[9], line_list[3], line_list[5], line_list[4])
        fastqseq_list.append(fastqseq)

        line = bamfile.readline()
        line_list = line.split('\t')

    for fastqseq in fastqseq_list:  # fastqseq_list存放着处理过的fastq序列
        for i in range(0, len(fastqseq)):  # 更新nt_list
            if fastqseq[i] == 'A':
                nt_list[i, 0] += 1
            elif fastqseq[i] == 'T':
                nt_list[i, 1] += 1
            elif fastqseq[i] == 'C':
                nt_list[i, 2] += 1
            elif fastqseq[i] == 'G':
                nt_list[i, 3] += 1
    nt_list=nt_list.astype(int)
    gene_name = gene_header.split('|')[1].replace('*', '_')
    if sum(sum(nt_list)) == 0:
        pass
    else:
       np.savetxt('nucleotideCounts/'+gene_name + '.csv', nt_list, fmt='%s',delimiter=';')  # 每个基因创建一个csv文件来记录各个位点的碱基数，因为csv文件的命名不能出现*|所以用了一些字符串操作函数


CDR_region = find_CDR()
csvfile_list = file_list('nucleotideCounts/', ext='.csv')
for csvfile in csvfile_list:
    freqs, depth = getAltFreq(csvfile[:-4])
    Mutcounts = freqs * np.transpose([depth, depth, depth, depth])
    Mutcounts = Mutcounts.astype(int)
    np.savetxt('MutCounts/'+ csvfile[:-4] + '.txt', Mutcounts, fmt='%s')
    plotMutations(freqs, csvfile[:-4], CDR_region[csvfile[:-4]])

referencefile.close()
bamfile.close()