

import os
from typing import SupportsComplex
import numpy as np
import pandas as pd
# import umap
import sys
import matplotlib.pyplot as plt

def removeduplicate(inputfilepath,outputfilpath):

    genenamelist=[]
    with open(outputfilpath,"w") as outputfile:
        with open(inputfilepath,"r") as inputfile:
            for line in inputfile.readlines():
                linedata = line.strip().split()
                genename = linedata[0]
                print(genename)
                if not genename in genenamelist:
                    genenamelist.append(genename)
                    outputfile.write(line)




def loadsamplelist(filedir):

    commongenelist=[]
    samplelist={}
    g = os.walk(filedir)
    for path,dir_list,file_list in g:  
        for file_name in file_list:  
            print(file_name)
            with open(os.path.join(path,file_name),"r") as matrixfile:
                filegenelist=[]
                isheader = True
                filesamleplist={}
                headerorder=[]
                for line in matrixfile.readlines():
                    if isheader:
                        linedata = line.strip().split("\t")
                        del(linedata[0])
                        for cellname in linedata:
                            filesamleplist[cellname]={}
                            headerorder.append(cellname)
                        isheader=False
                    else:
                        linedata = line.strip().split("\t")
                        genename = linedata[0].upper()
                        filegenelist.append(genename)
                        del(linedata[0])
                        for i in range(len(linedata)):
                            cellname = headerorder[i]
                            filesamleplist[cellname][genename]=float(linedata[i])
                if len(commongenelist)==0:
                    commongenelist=filegenelist
                else:
                    commongenelist = [i for i in commongenelist if i in filegenelist]                
                for cellname,genelist in filesamleplist.items():
                    samplelist[cellname] = genelist
                # free some memory 
                del(filesamleplist)

            print(len(samplelist))
    return {
        "commongenelist":commongenelist,
        "samplelist":samplelist
    }


def writeappendfile(fileoutput,fileoutputnoheader,samplelist):
    print("writing>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    with open(fileoutput,"w") as output:
        header = "GENE"
        for cellname in samplelist["samplelist"].keys():
            header+="\t"+cellname
        output.write(header+"\n")

        for genename in samplelist["commongenelist"]:
            # print(genename)
            pline=genename
            for cellname,genelist in samplelist["samplelist"].items():
                pline+="\t"+str(genelist[genename])
            output.write(pline+"\n")
    
    # with open(fileoutputnoheader,"w") as output:
    #     # header = "GENE"
    #     # for cellname in samplelist["samplelist"].keys():
    #     #     header+="\t"+cellname
    #     # output.write(header+"\n")

    #     for genename in samplelist["commongenelist"]:
    #         print(genename)
    #         # pline=genename
    #         pline =""
    #         for cellname,genelist in samplelist["samplelist"].items():
    #             pline+=str(genelist[genename])+"\t"
    #         output.write(pline[:-1]+"\n")

    

# def umapDR(filepath,outputpath):

#     dataset = pd.read_csv(filepath, delimiter = "\t",index_col=0)
#     print(dataset.shape)
#     # sys.exit()
#     reducer = umap.UMAP()

#     afterDRdata = reducer.fit_transform(dataset.T)

#     print(afterDRdata.shape)

#     # plt.scatter(
#     #     afterDRdata[:, 0],
#     #     afterDRdata[:, 1],
#     # )
    
#     # plt.title('UMAP Nash and Chow liver', fontsize=24)
#     # plt.show()

#     with open(outputpath,"w") as output:
#         for i in range(len(dataset.columns)):
#             pline=dataset.columns[i]+"\t"+str(afterDRdata[i,0])+"\t"+str(afterDRdata[i,1])+"\n"
#             output.write(pline)


#     # umap.plot.points(afterDRdata)

def identifiersingleCellSeurat(seuratmarkerfilepath,markergenfile,outputfile):
    makergenecelltypelist={}
    with open(markergenfile,"r") as hallmarkgene:
        for line in hallmarkgene.readlines():            
            linedata = line.strip().split("\t")
            genename = linedata[0].upper()
            celltype = linedata[1]
            if not celltype in makergenecelltypelist.keys():
                makergenecelltypelist[celltype] = [genename]
            else:
                makergenecelltypelist[celltype].append(genename) 

    seurateclusterlist = {}
    with open(seuratmarkerfilepath,"r") as seuratmarker:
        isheader = True
        for line in seuratmarker.readlines():
            if isheader:
                isheader = False
            else:
                linedata = line.replace("\"","").strip().split("\t")
                clusternumber = linedata[5]
                gene = linedata[6]
                if not clusternumber in seurateclusterlist.keys():
                    seurateclusterlist[clusternumber] = [gene]
                else:
                    seurateclusterlist[clusternumber].append(gene)
    with open(outputfile,"w") as output:
        # analysis
        for clusternumber,genelist in seurateclusterlist.items():
            for celltype,markgenelist in makergenecelltypelist.items():
                contains = 0
                for genename in markgenelist:
                    if genename in genelist:
                        contains+=1
                print(clusternumber+"\t"+celltype+"\t"+str(contains))
                output.write(clusternumber+"\t"+celltype+"\t"+str(contains)+"\n")

           


def identifiersingleCell(filepath,hallmarkgenefile,outputfilepath,threshold=1):

    makergenecelltypelist={}
    celltypelist = {}
    index = 1
    celltypelist["nan"] = 0

    with open(hallmarkgenefile,"r") as hallmarkgene:
        for line in hallmarkgene.readlines():            
            linedata = line.strip().split("\t")
            genename = linedata[0].upper()
            celltype = linedata[1]
            if not celltype in makergenecelltypelist.keys():
                makergenecelltypelist[celltype] = [genename]
                celltypelist[celltype]=index
                index+=1
            else:
                makergenecelltypelist[celltype].append(genename)

    sampletype = {}
    samplegenelist = {}
    index = 0

    filesamleplist={}
    with open(filepath,"r") as matrixfile:
        isheader = True        
        headerorder=[]
        for line in matrixfile.readlines():
            if isheader:
                linedata = line.strip().split("\t")
                del(linedata[0])
                for cellname in linedata:
                    filesamleplist[cellname]=[]
                    headerorder.append(cellname)
                isheader=False
            else:
                linedata = line.strip().split("\t")
                genename = linedata[0].upper()
                samplegenelist[genename]=index
                index+=1
                del(linedata[0])
                for i in range(len(linedata)):
                    cellname = headerorder[i]
                    filesamleplist[cellname].append(float(linedata[i]))

    for samplename,samplegenevaluelist in filesamleplist.items():
        for celltype,markergenelist in makergenecelltypelist.items():
            istype = True
            for genename in markergenelist:
                if genename not in samplegenelist or samplegenevaluelist[samplegenelist[genename]]<threshold:
                    istype=False
                    break
            if istype:
                sampletype[samplename]=celltype
                break
        if samplename not in sampletype.keys():
            sampletype[samplename] = "nan"
        print(samplename+": "+sampletype[samplename])
    
    with open(outputfilepath,"w") as outputfile:
        for samplename,celltype in sampletype.items():
            outputfile.write(samplename+"\t"+celltype+"\t"+str(celltypelist[celltype])+"\n")

def selectsubgroupofMARK(originalfilepath,markergenfile,outputfile):


    makergenecelltypelist=[]
  

    with open(markergenfile,"r") as hallmarkgene:
        for line in hallmarkgene.readlines():            
            linedata = line.strip().split("\t")
            genename = linedata[0].upper()

            if not genename in makergenecelltypelist:
                makergenecelltypelist.append(genename)
    

    with open(outputfile,"w") as output:
        with open(originalfilepath,"r") as originalfile:
            isheader = True        
            for line in originalfile.readlines():
                if isheader:
                    output.write(line)
                    isheader=False
                else:
                    linedata = line.strip().split("\t")
                    genename = linedata[0].upper()
                    if genename in makergenecelltypelist:
                        output.write(line)


                    

# if __name__=="__main__":
    
    # filedir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/SingleCell/chow/"

    # # fileoutput = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/test.txt"
    # fileoutput = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/combined.txt"
    # outputfilpath =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/combined_removedup.txt"
     
    # hallmarkgenefile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/makertest.txt"

    # umapDRoutput = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/combinedUMAPdr.txt"

    # outputfilepath =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/label.txt"

    # samplelist = loadsamplelist(filedir)

    # writeappendfile(fileoutput,fileoutputnoheader,samplelist)

    # umapDR(fileoutput,umapDRoutput)
 
    # identifiersingleCell(fileoutput,hallmarkgenefile,outputfilepath,threshold=1)

    # removeduplicate(fileoutput,outputfilpath)
 
    # seuratmarkerfilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/3-9/seurat_marker_1.txt"
    # markergenfile="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/makertest.txt"

    # outputfilepath =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/3-9/label_1.txt"

    # identifiersingleCellSeurat(seuratmarkerfilepath,markergenfile,outputfilepath)
                        
    # outputfile="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE119340/combined_subgroup.txt"
    # selectsubgroupofMARK(fileoutput,hallmarkgenefile,outputfile)



