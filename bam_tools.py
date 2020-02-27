import re
import numpy as np

def align(reference_seq, seq, pos, CIGAR, quality):
    CIGAR_list = re.findall(r'[0-9]+|[A-Z]+', CIGAR)
    pos = int(pos)
    seq_new = ''
    quality_new = ''
    for i in range(0, len(CIGAR_list) // 2):
        if CIGAR_list[i*2+1] == 'S':
            seq = seq[int(CIGAR_list[2*i]):]
            quality = quality[int(CIGAR_list[2*i]):]
        elif CIGAR_list[2*i+1] == 'M':
            seq_new = seq_new+seq[0:int(CIGAR_list[2*i])]
            quality_new = quality_new + quality[0:int(CIGAR_list[2*i])]
            seq = seq[int(CIGAR_list[2*i]):]
            quality = quality[int(CIGAR_list[2*i]):]
        elif CIGAR_list[2*i+1] == 'D':
            seq_new = seq_new+'D'
            quality_new = quality_new + 'D'
        elif CIGAR_list[2*i+1] == 'I':
            seq = seq[int(CIGAR_list[2*i]):]
            quality = quality[int(CIGAR_list[2*i]):]
        elif CIGAR_list[len(CIGAR_list)-1] == 'S':
            pass
    seq_new = (pos-1)*'N'+seq_new
    quality_new = (pos-1)*'N'+quality_new
    return (seq_new, quality_new)

# print(align(reference_seq, seq, pos, CIGAR, quality))


def find_mismatch_index(seqs):
    mismatch_index = []
    seq1_list = list(seqs[0])
    seq2_list = list(seqs[1])
    for i in range(0, len(seq2_list)):
        if seq1_list[i] != seq2_list[i] and seq2_list[i] != 'N':
            mismatch_index.append(i)
    return mismatch_index


def getAltFreq(sample):
    # 这个函数可以统计各个位点的突变的碱基频率，输入的是文件的名字（文件在nucleotideCounts文件夹下），
    # 输出的是频率（n行4列的数组）与深度（一位数组）
    csv = 'nucleotideCounts/' + sample + '.csv'
    Counts = []
    with open(csv, 'r') as f:
        for line in f:
            Counts.append([int(c) for c in line.rstrip(';\n').split(';')])

    # Compute nucleotide frequencies and depth
    Depths = []
    Freqs = []
    for count in Counts:
        Depths.append(sum(count))
        Freqs.append([c / sum(count) if sum(count)>0 else 1 for c in count])
    Freqs = np.array(Freqs)

    # Filter reference nucleotide at each position
    (n, m) = np.shape(Freqs)
    for i in range(n):
        for j in range(m):
            if Freqs[i, j] > 0.5:
                Freqs[i, j] = 0
    return (Freqs, Depths)

def find_CDR(filename="Human IGHV F+ORF+in-frame P含..txt"):
    referencefile = open(filename)
    reference_genelist = referencefile.readlines()
    CDR_dict = {}
    for i in range(0, len(reference_genelist) // 2):  # 2i是基因名，2i+1是基因
        CDR1_start = 27
        CDR1_end = 113
        CDR2_start = 167
        CDR2_end = 194
        fastqseq_list = []
        gene_header = reference_genelist[2 * i].strip()
        reference_gene = reference_genelist[2 * i + 1].strip()
        gene_name = gene_header.split('|')[1].replace('*', '_')
        for i in range(len(reference_gene)):
            if reference_gene[i] == '.' and i < CDR1_start:
                CDR1_start = CDR1_start - 1
                CDR1_end = CDR1_end - 1
                CDR2_start = CDR2_start - 1
                CDR2_end = CDR2_end - 1
            elif reference_gene[i] == '.' and i >= CDR1_start and i <= CDR1_end:
                CDR1_end = CDR1_end - 1
                CDR2_start = CDR2_start - 1
                CDR2_end = CDR2_end - 1
            elif reference_gene[i] == '.' and i < CDR2_start:
                CDR2_start = CDR2_start - 1
                CDR2_end = CDR2_end - 1
            elif reference_gene[i] == '.' and i >= CDR2_start and i <= CDR2_end:
                CDR2_end = CDR2_end - 1
        CDR_dict[gene_name] = [CDR1_start,CDR1_end,CDR2_start,CDR2_end]
    return CDR_dict


"""gene='CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA'
freqs, depth = getAltFreq('IGHV1-18_01')
freqs_mutant = []
for freq in freqs:
    freqs_mutant.append(sum(freq))
mutant_index = np.nonzero(freqs_mutant)
print(mutant_index)

Mutcounts = freqs * np.transpose([depth, depth, depth, depth])
print(freqs[range(56,75)])
CDR_region = find_CDR()
print(CDR_region['IGHV1-18_01'])"""