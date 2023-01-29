# This program trying to make the searching space smaller.
# By calcualte the pure degree and choose top x coreelation link  
    







from scipy import stats
import itertools
from multiprocessing import Process
import operator


def calculatecorrelation_mult(correlationType,xlist,ylist):
    returndict={}
    
    if correlationType == "pearsonr":
        cor,pval = stats.pearsonr(xlist,ylist)
        returndict['cor'] = cor
        returndict['pval'] = pval
    elif correlationType == "spearmanr":
        cor,pval =stats.spearmanr(xlist,ylist)
        returndict['cor'] = cor
        returndict['pval'] = pval
    elif correlationType == "kendalltau":
        cor,pval = stats.kendalltau(xlist,ylist)
        returndict['cor'] = cor
        returndict['pval'] = pval
    
    return returndict

def filteroutLinks_Correlation(ExpressionGeneCount,linkslist, toplinkNumber = 5 ,correlationType = "pearsonr",outputfile=False):

    totallink={}
    
    for link in linkslist:
        gene1 = link["gene1"]
        gene2 = link["gene2"]
        xlist = ExpressionGeneCount[gene1]
        ylist = ExpressionGeneCount[gene2]
        rlist = calculatecorrelation_mult(correlationType,xlist,ylist)
        if gene1 not in totallink.keys():
            totallink[gene1]={
               gene2:abs(rlist["cor"]), 
            }
        else:
            totallink[gene1][gene2]=abs(rlist["cor"])
        
        if gene2 not in totallink.keys():
            totallink[gene2]={
               gene1:abs(rlist["cor"]), 
            }
        else:
            totallink[gene2][gene1]=abs(rlist["cor"])

    newlinkslist =[]

    for gene1, endlist in totallink.items():
        newA = dict(sorted(endlist.items(), key=operator.itemgetter(1), reverse=True)[:toplinkNumber])
        for gene2,value in newA.items():
            newlink={
                    "gene1":gene1,
                    "gene2":gene2
                }
            newlinkslist.append(newlink)
    print(len(newlinkslist))
    if  outputfile:
        with open(outputfile,"w") as output:
            for link in newlinkslist:
                pline = link["gene1"]+"\t"+link["gene2"]+"\n"
                output.write(pline)
    return newlinkslist


def CohortGenepureDegree(GeneCount):
    ConsistenceData = {}
           
    for genename,Valuelist in GeneCount.items():
        positive=0
        negative =0
        for value in Valuelist:
            if value >0:
                positive+=1
            if value <0:
                negative+=1
        Log2ConsistDegree=max([positive,negative])/len(Valuelist)
        ConsistenceData[genename]=Log2ConsistDegree
       
    return ConsistenceData







def filteroutLinks_log2(Log2GeneCount,linkslist,puredegree = 0.7,outputfile=False):
    # filterGene
    GenePure = CohortGenepureDegree(Log2GeneCount)

    # star filter
    newlinkslist=[]
    
    for link in linkslist:
        gene1 = link["gene1"]
        gene2 = link["gene2"]
        if gene1 in GenePure.keys() and  gene2 in GenePure.keys():
            if GenePure[gene1]>= puredegree and  GenePure[gene2]>= puredegree :
                newlink={
                    "gene1":gene1,
                    "gene2":gene2
                }
                newlinkslist.append(newlink)
    print(len(newlinkslist))
    if  outputfile:
        with open(outputfile,"w") as output:
            for link in newlinkslist:
                pline = link["gene1"]+"\t"+link["gene2"]+"\t"+str(GenePure[link["gene1"]])+"\t"+str(GenePure[link["gene2"]])+"\n"
                output.write(pline)
    return newlinkslist



def generateGraph(linkfilepath):

    genepairs = []

    with open(linkfilepath,"r") as genepairsfile:
        for line in genepairsfile.readlines():
            linkarray = {}
            linedata= line.strip().split("\t")
            linkarray['gene1'] = linedata[0].upper()
            linkarray['gene2'] = linedata[1].upper()
            genepairs.append(linkarray)


    graph ={}
    returngraph = {}

    for genepair in genepairs:
        if genepair["gene1"] not in graph.keys():
            graph[genepair["gene1"]]=[]
            graph[genepair["gene1"]].append(genepair["gene2"])
        else:
            graph[genepair["gene1"]].append(genepair["gene2"])

        if genepair["gene2"] not in graph.keys():
            graph[genepair["gene2"]]=[]
            graph[genepair["gene2"]].append(genepair["gene1"])
        else:
            graph[genepair["gene2"]].append(genepair["gene1"])
    
    for key,linkednodelist in graph.items():
        returngraph[key]=set(linkednodelist)
    
    return returngraph

def generateRegGraph(linkfilepath):

    genepairs = []

    with open(linkfilepath,"r") as genepairsfile:
        for line in genepairsfile.readlines():
            linkarray = {}
            linedata= line.strip().split("\t")
            linkarray['gene1'] = linedata[0].upper()
            linkarray['gene2'] = linedata[2].upper()
            genepairs.append(linkarray)
    
    graph ={}
    returngraph = {}

    for genepair in genepairs:
        if genepair["gene1"] not in graph.keys():
            graph[genepair["gene1"]]=[]
            graph[genepair["gene1"]].append(genepair["gene2"])
        else:
            graph[genepair["gene1"]].append(genepair["gene2"])
            
    
    for key,linkednodelist in graph.items():
        returngraph[key]=set(linkednodelist)
    
    return returngraph






if __name__ == "__main__":

    correlationType = "pearsonr"
    puredegree=0.7
    toplinkNumber=2

    geneRatiopath = "./Cancer_data/TCGA_COAD_log2ratio.txt"
    geneExpressionpath = "./Cancer_data/TCGA_COAD_raw.txt"
    linkfilepath  = "DBlinks.txt"

    # read ratio
    GeneCount_log2={}
    # read routeScore table
    with open (geneRatiopath,"r") as geneLogfile:
        title = True
        for line in geneLogfile.readlines():
            if title:
                title =False
            else:
                linedata = line.strip().split("\t")
                genename = linedata[0].upper()
                del(linedata[0])
                valuecontent = [float(i) for i in linedata]
                GeneCount_log2[genename] = valuecontent
    print("read ratio")

    # read raw
    GeneCount_raw={}
    # read routeScore table
    with open (geneExpressionpath,"r") as generawfile:
        title = True
        for line in generawfile.readlines():
            if title:
                title =False
            else:
                linedata = line.strip().split("\t")
                genename = linedata[0].upper()
                del(linedata[0])
                valuecontent = [float(i) for i in linedata]
                GeneCount_raw[genename] = valuecontent
    print("read raw")
    #readlinks
    linkslist=[]
    with open (linkfilepath,"r") as linkfile:
        for line in linkfile.readlines():
                linedata = line.strip().split("\t")
                newlink={
                    "gene1":linedata[0],
                    "gene2":linedata[1]
                }
                linkslist.append(newlink)
    print(len(linkslist))


    logfiterlinks = filteroutLinks_log2(GeneCount_log2,linkslist,puredegree = 0.7)
    corfiterlinks = filteroutLinks_Correlation(GeneCount_raw,logfiterlinks, toplinkNumber,correlationType = "pearsonr",outputfile="DBlinks_"+str(puredegree)+"_"+str(toplinkNumber)+".txt")





#     Filterlinkfilepath  = "DBlinks_log2_filter_"+str(puredegree)+".txt"

#     # # calculate gene corhort correlation
#     # GeneCount={}
#     # # read routeScore table
#     # with open (geneRatiopath,"r") as geneLogfile:
#     #     title = True
#     #     for line in geneLogfile.readlines():
#     #         if title:
#     #             title =False
#     #         else:
#     #             linedata = line.strip().split("\t")
#     #             genename = linedata[0].upper()
#     #             del(linedata[0])
#     #             valuecontent = [float(i) for i in linedata]
#     #             GeneCount[genename] = valuecontent

#     # genepairs = []          
#     # print("Write genepairsfile_cor")
#     # with open(linkfilepath,"r") as genepairsfile:
#     #     with open(Filterlinkfilepath,"w") as Filterlinkfile:
#     #         i =0
#     #         q=0
            
#     #         lines = genepairsfile.readlines()
#     #         for line in lines:
#     #             i+=1
#     #             prob = i/len(lines)*100
                                               
#     #             linedata= line.strip().split("\t")
#     #             gene1 = linedata[0].upper()
#     #             gene2 = linedata[1].upper()
#     #             if gene1 in GeneCount.keys() and  gene2 in GeneCount.keys():
#     #                 rlist = calculatecorrelation_mult(correlationType, GeneCount[gene1], GeneCount[gene2])
#     #                 if abs(rlist['cor']) >= linkcorrelation and rlist['pval'] <= linkPval:
#     #                     printline = gene1+"\t"+gene2+"\t"+str(round(rlist['cor'],3))+"\t"+str(round(rlist['pval'],3))+"\n"
#     #                     print(prob)
#     #                     q+=1

#     #                     Filterlinkfile.write(printline)

#     # print("orinal:",str(len(lines)))
#     # print("New:",str(q))

    
#     GeneCount={}
#     # read routeScore table
#     with open (geneRatiopath,"r") as geneLogfile:
#         title = True
#         for line in geneLogfile.readlines():
#             if title:
#                 title =False
#             else:
#                 linedata = line.strip().split("\t")
#                 genename = linedata[0].upper()
#                 del(linedata[0])
#                 valuecontent = [float(i) for i in linedata]
#                 GeneCount[genename] = valuecontent
#     filterroutegraph_log2(GeneCount,Filterlinkfilepath,linkfilepath,puredegree)      




