import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from bam_tools import *


def file_list(dirname, ext='.csv'):
    """获取目录下所有特定后缀的文件
    @param dirname: str 目录的完整路径
    @param ext: str 后缀名, 以点号开头
    @return: list(str) 所有子文件名(不包含路径)组成的列表
    """
    return list(filter(lambda filename: os.path.splitext(filename)[1] == ext,os.listdir(dirname)))


def plotMutations(Freqs, figName, display=False, region=0):
    ### This function returns NumPy arrays of prediction intervals.
    ### Parameters:    Test: Numpy array [length(gene) x 4]
    ###                        Alternative nucleotide frequencies in sample as returned by getAltFreq
    ###                Testd: Numpy array [length(gene) x 1]
    ###                        Position-wise depth in sample as returned by getAltFreq
    ###                gene: string
    ###                        File "gene.fasta" must be in nucleotideCounts directory
    ###                region: range
    ###                        Specify the region of gene sequence.
    ###                figName: string
    ###                        File name to save generated figures

    # Mutcounts = Freqs * np.transpose([Depths, Depths, Depths, Depths])
    colors = ['#57bddb', '#ea4848', '#3bb971', '#f39c12']

    # Plot position-wise mismatch frequency
    plt.plot(Freqs[:, 0] + Freqs[:, 1] + Freqs[:, 2] + Freqs[:, 3], 'b')
    plt.rcParams.update({'font.size': 10})
    plt.axis([0, len(Freqs), 0, 0.2])
    plt.xlabel('Position (bp)')
    plt.ylabel('Mutation frequency')
    if display:
        plt.show()
    plt.savefig('MutCounts/'+figName + '.pdf', bbox_inches='tight')
    plt.close()

    # Plot position-wise mismatch frequency for each nucleotide
    ymax = 0.2  # 突变频率上限设定为0.2
    f, axarr = plt.subplots(2, 2)
    plt.rcParams.update({'font.size': 10})
    axarr[0, 0].plot(Freqs[:, 0], color=colors[0])  # 第0列是A
    axarr[0, 0].set_ylim([0, ymax])
    axarr[0, 0].set_xlim([0, len(Freqs)])
    axarr[0, 0].set_title('Mut. to A',color=colors[0])
    axarr[0, 1].plot(Freqs[:, 1], color=colors[1])  # 第1列是T
    axarr[0, 1].set_ylim([0, ymax])
    axarr[0, 1].set_xlim([0, len(Freqs)])
    axarr[0, 1].set_title('Mut. to T',color=colors[1])
    axarr[1, 0].plot(Freqs[:, 2], color=colors[2])  # 第2列是C
    axarr[1, 0].set_ylim([0, ymax])
    axarr[1, 0].set_xlim([0, len(Freqs)])
    axarr[1, 0].set_title('Mut. to C',color=colors[2])
    axarr[1, 1].plot(Freqs[:, 3], color=colors[3])  # 第3列是G
    axarr[1, 1].set_ylim([0, ymax])
    axarr[1, 1].set_xlim([0, len(Freqs)])
    axarr[1, 1].set_title('Mut. to G',color=colors[3])
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    axarr[1, 1].set_xlabel('Position (bp)')
    axarr[1, 0].set_xlabel('Position (bp)')
    axarr[1, 0].set_ylabel('Mut. frequency')
    axarr[0, 0].set_ylabel('Mut. frequency')
    if display:
        plt.show()
    plt.savefig('MutCounts/'+figName + '_bases.pdf', bbox_inches='tight')
    plt.close()
    return


