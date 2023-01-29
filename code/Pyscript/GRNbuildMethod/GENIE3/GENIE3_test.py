from GENIE3 import *
import sklearn
import numpy as np


import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

def normalization(data):
    _range = np.max(data) - np.min(data)
    return (data - np.min(data)) / _range



if __name__ == "__main__":
    index = 9
    # for index in range(8):

    data = np.loadtxt("C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/randomselect/test4/"+str(index)+"_test.txt")
    data = data.transpose()
    data = normalization(data)
    f = open("C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/randomselect/test4/"+str(index)+"_genename.txt")
    # gene_names = f.readline()
    # f.close()

    gene_names=[]
    for line in f.readlines():
        gene_names.append(line.strip()) 
    f.close()
    # gene_names = gene_names.strip('\n').split('\t')

    VIM = GENIE3(data,nthreads=16)
    # # get_link_list(VIM)
    # # get_link_list(VIM,gene_names=gene_names)

    get_link_list(VIM,gene_names=gene_names,file_name="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/randomselect/test4/result/"+str(index)+"_result.txt")

    # ex_matrix = pd.read_csv('C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/GRNbuildMethod/GENIE3/STAGE 0vsIIIC_5_5.txt', sep='\t')


    # network = grnboost2(expression_data=ex_matrix)

    # network.to_csv('C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/GRNbuildMethod/GENIE3/output.tsv', sep='\t', index=False, header=False)




