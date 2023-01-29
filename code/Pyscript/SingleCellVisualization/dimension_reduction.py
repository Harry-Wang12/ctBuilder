## importing the required packages
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas
from sklearn import datasets
from sklearn import decomposition
from sklearn import manifold
import umap
from time import time
import random
import os, sys

'''
@author Pujan Joshi
@date Feb 01, 2019
This code runs PCA, UMAP and tSNE.
This code can also draw scatter plot from results of PCA, UMAP and tSNE.
This code also writes resuts (X-Y info) to the files.

****library install scripts***

python -m pip install -U pip
python -m pip install -U matplotlib
python -m pip install pandas
python -m pip install -U scikit-learn
python -m pip install umap-learn

****library install scripts end here***

'''
work_dir = ""
inputFilename = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/GSE155182_Normal_7m_liver.txt' # GSE115469_P5.csv
outputFilename = ''

RUN_PCA =  True #True|False
RUN_UMAP = True #True|False
RUN_TSNE = True #True|False
SHOW_PLOT = True #True|False
TRANSPOSE = False #True|False #Transpose is needed if samples are rows
n_iter = 700
delimiter = '\t'; #'\t' or ','

#####
##DO NOT CHANGE ANYTHING BELOW THIS LINE
#####

if work_dir != "":
    os.chdir(work_dir)
print("Current Working Directory " , os.getcwd())

help_msg = '''command line usage: 
python dimension_reduction.py 
-i <input file full path> (required)
-o <output file full path> (default '', will auto-generate)
-niter 700  (default)
-delimiter \\t (default) || , 
-runPCA True || False (default)
-runUMAP True || False (default)
-runTSNE True || False (default)
-showPlot True || False (default)
-transpose True || False (default)
'''
def parse_sys_argv(keys, help_msg=help_msg):
    args = sys.argv
    arguments = {}
    for key in keys:
        if key in args:
            arguments[key] = args[args.index(key) + 1]
    if len(arguments) == 0:
        print()
        print(help_msg)
        print()
        sys.exit()
    return arguments

# # args_keys = ["-i", "-o", '-niter', '-delimiter', '-runPCA', '-runUMAP', '-runTSNE', '-showPlot', '-transpose', '-addLabel', '-labelFile', '-labelIndex']
# # args = parse_sys_argv(args_keys)

# # if '-i' in args:
#     # inputFilename = args['-i']
# inputFilename = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/combined_subgroup.txt"
# # if '-o' in args:
# outputFilename = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/combined_subgroup_test.txt"
# # if '-niter' in args:
# n_iter = 200
# if '-delimiter' in args:
# delimiter = args['-delimiter']
# if '-runPCA' in args:
#     # runPCA = args['-runPCA']
#     # RUN_PCA = False
#     # if runPCA.upper() == 'TRUE':
#     RUN_PCA = True
# if '-runUMAP' in args:
#     # runUMAP = args['-runUMAP']
#     # RUN_UMAP = False
#     # if runUMAP.upper() == 'TRUE':
#     RUN_UMAP = True        
# if '-runTSNE' in args:
#     # runTSNE = args['-runTSNE']
#     # RUN_TSNE = False
#     # if runTSNE.upper() == 'TRUE':
#     RUN_TSNE = True        
# if '-showPlot' in args:
#     # showPlot = args['-showPlot']
#     # SHOW_PLOT = False
#     # if showPlot.upper() == 'TRUE':
#     SHOW_PLOT = True

# if '-transpose' in args:
#     transpose = args['-transpose']
#     TRANSPOSE = False
#     if transpose.upper() == 'TRUE':
#         TRANSPOSE = True        
print("start")
t = round(time())
rand_int = random.randint(1,100)
o_filename_suffix = "_" + str(t) + "_" + str(rand_int)

if outputFilename == '':
    outputFilename = inputFilename

dot_pos = outputFilename.rfind(".")
if dot_pos >= 0:
    outputFilename = outputFilename[: dot_pos]

o_filename = outputFilename + o_filename_suffix

## Function to Scale and visualize the embedding vectors
WRITE_TO_FILE = True

def plot_embedding(X, title=None):
    x_min, x_max = np.min(X, 0), np.max(X, 0)
    X = (X - x_min) / (x_max - x_min)     
    plt.figure()
    ax = plt.subplot(111)
    for i in range(X.shape[0]):
        plt.text(X[i, 0], X[i, 1], ".", 
                 fontdict={'weight': 'bold', 'size': 9})
    
    plt.xticks([]), plt.yticks([])
    if title is not None:
        plt.title(title)
if TRANSPOSE:
    with open(inputFilename) as file:
        lis = [x.replace('\n', '').split('\t') for x in file]

    data = np.array(lis).T
    df = pandas.DataFrame(data=data[1:,:], index=data[1:,0], columns=data[0,:])
else:
    df = pandas.read_csv(inputFilename, delimiter=delimiter)

cols = list(df.columns)[1:]
X = df.values
X = X[:, 1:]
X = np.transpose(X)
cols_t = np.array([cols])
if RUN_PCA:
    ## Computing PCA
    print()
    print("Computing PCA projection...")
    t0 = time()
    X_pca = decomposition.TruncatedSVD(n_components=2).fit_transform(X)
    #X_pca = decomposition.PCA(n_components=2).fit(X)
    tt_pca = round(time() - t0, 2)
    print("time taken by PCA:", tt_pca, "seconds")    
    if SHOW_PLOT:
        print("generating PCA plot...")
        plot_title = "PCA projection (time %.2fs)" % (tt_pca)
        plot_embedding(X_pca, plot_title)

    if WRITE_TO_FILE:
        pca_filename = o_filename + '_pca.txt'
        print("writing PCA results to file named", pca_filename)
        X_pca_2 = np.concatenate((cols_t.T, X_pca), axis=1)
        np.savetxt(pca_filename, X_pca_2, delimiter=delimiter, fmt='%s')

if RUN_UMAP:
    ## Computing UMAP
    print()
    print("Computing UMAP projection...")
    t0 = time()
    umap_dr = umap.UMAP()
    X_umap = umap_dr.fit_transform(X)
    tt_umap = round(time() - t0, 2)
    print("time taken by UMAP:", tt_umap, "seconds")
    if SHOW_PLOT:
        print("generating umap scatter plot...")
        plot_title = "UMAP Projection (time %.2fs)" % (tt_umap)
        plot_embedding(X_umap, plot_title)
    if WRITE_TO_FILE:
        umap_filename = o_filename + '_umap.txt'
        print("writing UMAP results to file named", umap_filename)
        X_umap_2 = np.concatenate((cols_t.T, X_umap), axis=1)
        np.savetxt(umap_filename, X_umap_2, delimiter=delimiter, fmt='%s')
        
if RUN_TSNE:
    ## Computing t-SNE
    print()
    print("Computing t-SNE with n_iter", n_iter, "..")
    tsne = manifold.TSNE(init='pca', random_state=1, n_iter=n_iter)
    #tsne = manifold.TSNE()
    t0 = time()
    print("Transformation..")
    X_tsne = tsne.fit_transform(X)
    tt_tsne = round(time() - t0, 2)
    print("time taken by tSNE:", tt_tsne, "seconds")
    if SHOW_PLOT:
        print("generating t-SNE plot...")
        plot_title = "t-SNE (n_iter %d)(time %.2fs)" % (n_iter, tt_tsne)
        plot_embedding(X_tsne, plot_title)
    if WRITE_TO_FILE:
        tsne_filename = o_filename + '_tsne_n' + str(n_iter) + '.txt'
        print("writing tSNE results to file named", tsne_filename)
        X_tsne_2 = np.concatenate((cols_t.T, X_tsne), axis=1)
        np.savetxt(tsne_filename, X_tsne_2, delimiter="\t", fmt='%s')

if SHOW_PLOT:
    print("displaing plots now ..")
    plt.show()

