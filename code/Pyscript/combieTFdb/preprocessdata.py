
from scipy.stats import ranksums
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.spatial.distance import cdist
import numpy as np
import matplotlib.pyplot as plt
import sys


def findthedifferentgene(filepath1,filepath2,degenefilepath,threshold):
    content1={}
    content2={}
    genelist = []
    degenelist=[]
    titleline = "Gene"
    # load file1

    with open(filepath1,'r') as file1:
        istitle=True

        for line in file1.readlines():
            linedata = line.strip().split("\t")
            if istitle:
                del(linedata[0])
                for titlename in linedata:
                    titleline+="\t"+titlename
                istitle=False
            else:
                
                genename = linedata[0].upper()
                del(linedata[0])
                content1[genename]=[]
                for genevalue in linedata:
                    content1[genename].append(float(genevalue))
                if not genename in genelist:
                    genelist.append(genename)
    
    with open(filepath2,'r') as file2:
        istitle=True

        for line in file2.readlines():
            linedata = line.strip().split("\t")
            if istitle:
                del(linedata[0])
                for titlename in linedata:
                    titleline+="\t"+titlename
                istitle=False
            else:
                
                genename = linedata[0].upper()
                del(linedata[0])
                content2[genename]=[]
                for genevalue in linedata:
                    content2[genename].append(float(genevalue))
                if not genename in genelist:
                    genelist.append(genename)
    with open(degenefilepath,'w') as outputfile:
        for genename in genelist:
            # print(genename)
            if genename in content1.keys() and genename in content2.keys():
                staticvalue, pvalue  = ranksums(content1[genename], content2[genename])
                if pvalue<= threshold:
                    # print(genename)
                    outputfile.write(genename+'\n')
                    degenelist.append(genename)
    
    return degenelist

    # with open(outputfile1path,'w') as outputfile:
    #     outputfile.write(titleline+'\n')
    #     for genename in genelist:
    #         if genename in content1.keys() and genename in content2.keys():
    #             staticvalue, pvalue  = ranksums(content1[genename], content2[genename])
    #             if pvalue<= threshold:
    #                 print(genename,pvalue)
    #                 pline=genename
    #                 for genevalue in content1[genename]:
    #                     pline+="\t"+str(genevalue)
    #                 # for genevalue in content2[genename]:
    #                 #     pline+="\t"+str(genevalue)
    #                 outputfile.write(pline+'\n')

    # with open(outputfile2path,'w') as outputfile:
    #     outputfile.write(titleline+'\n')
    #     for genename in genelist:
    #         if genename in content1.keys() and genename in content2.keys():
    #             staticvalue, pvalue  = ranksums(content1[genename], content2[genename])
    #             if pvalue<= threshold:
    #                 print(genename,pvalue)
    #                 pline=genename
    #                 # for genevalue in content1[genename]:
    #                 #     pline+="\t"+str(genevalue)
    #                 for genevalue in content2[genename]:
    #                     pline+="\t"+str(genevalue)
    #                 outputfile.write(pline+'\n')


def kmeanselbowopt(kmeanslist,zeroline,maxclustering,minchange):

    returnlist = []
    X = np.array(kmeanslist)

    lastdist = sys.maxsize
    valuelabel=0
    centervalues=0
    

    for k in range(1,maxclustering):
        kmeanModel = KMeans(n_clusters=k).fit(X)
        kmeanModel.fit(X)
        nowdist = sum(np.min(cdist(X, kmeanModel.cluster_centers_,'euclidean'), axis=1)) / X.shape[0]
        valuelabel = list(kmeanModel.labels_)
        centervalues = list(kmeanModel.cluster_centers_)
        if lastdist - nowdist <= minchange and nowdist<=1:
            print(k,nowdist)
            
            break
        else:
            lastdist = nowdist
    
    
    index = 0
    
    for value in zeroline:
        if value:
            returnlist.append(0)
        else:
            returnlist.append(int(centervalues[valuelabel[index]]))
            index+=1
    
    return returnlist


def kmeanselbowplot(inputfilepath,outputfilepath,maxclustering=50,minchange = 0.1 ):
    
    with open (outputfilepath,'w') as outputfile:
        with open(inputfilepath,'r') as inputfile:
            istitle = True
            for line in inputfile.readlines():
                if istitle:
                    outputfile.write(line)
                    istitle = False
                else:
                    linedata = line.strip().split("\t")
                    genename = linedata[0]
                    pline = genename
                    del(linedata[0])
                    kmeansline = []
                    zeroline = []
                    uniquelist = []
                    for value in linedata:
                        if value =="0":
                            zeroline.append(True)
                        else:
                            zeroline.append(False)
                            if value not in uniquelist:
                                uniquelist.append(float(value))                                
                            kmeansline.append([float(value)])                    
                    if len(kmeansline)>maxclustering: 
                        if max(uniquelist)-min(uniquelist)>5:
                            print(genename)
                            for kmeansvalue in kmeanselbowopt(kmeansline,zeroline,min([maxclustering,len(uniquelist)]),minchange):
                                pline+="\t"+str(kmeansvalue)
                            pline+='\n'
                        else:
                            pline=line                        
                        outputfile.write(pline)






if __name__=="__main__":

    randomsamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/singlecellsimulation/counts_0.txt"

    kmeansrandomsamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/singlecellsimulation/kmeans_counts_0.txt"


    kmeanselbowplot(randomsamples,kmeansrandomsamples )

    # X=np.array([[1],[2],[2],[3],[4],[5],[5123],[123],[4],[5],[6],[3],[12],[36],[123],[14],[512],[6],[123],[623],[6],[1],[123],[53214],[7],[7],[8]])
    # k = 8
    # kmeanModel = KMeans(n_clusters=k).fit(X)
    # kmeanModel.fit(X)
    # print(kmeanModel.labels_)
    # print(kmeanModel.cluster_centers_)

    # print(sum(np.min(cdist(X, kmeanModel.cluster_centers_,'euclidean'), axis=1)) / X.shape[0])

#     filepath1="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/AD5_AD6/OPC.txt"
#     filepath2="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/Ct1_Ct2/OPC.txt"
#     threshold=0.1
#     outputfilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/combinedranksum/0.1/GSE138852_AD_control_OPC_small_"+str(threshold)+".txt"
#     findthedifferentgene(filepath1,filepath2,outputfilepath,threshold)



