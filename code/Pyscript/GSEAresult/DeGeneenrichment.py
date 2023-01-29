

from statistics import mean
import os
import sys


def loadroutes(routefile):
    paths = {}
    with open(routefile,'r') as routeinput:
        for line in routeinput.readlines():
            linedata = line.strip().split("\t")
            pathname = linedata[0]
            del(linedata[0])
            del(linedata[0])
            if len(linedata)>2:
                paths[pathname] = linedata
    
    return paths


def loadratiofile(filepath):
    CPlist={}
    with open(filepath,'r') as inputfile:
        istitle = True
        for line in inputfile.readlines():
            if not line.startswith("!"):
                if istitle:
                    istitle=False
                else:
                    linedata=line.strip().split('\t')
                    
                    genename=linedata[0]
                    ratio  = float(linedata[1])
                    # controlValue = float(linedata[2])
                    # testValue = float(linedata[3])
                    # controlZ = float(linedata[5])
                    # testZ = float(linedata[7])
                    CPlist[genename]={
                        'ratio':ratio,
                        # 'controlValue':controlValue,
                        # 'testValue':testValue,
                    # 'controlZ':controlZ,
                    # 'testZ':testZ                   

                    }
    return CPlist

def orderDEgene(CPlist,TC=10,TT=10):
    rawCPlist = {}
    for genename,valueinformation in CPlist.items():
        # if valueinformation['controlValue'] >=TC and valueinformation['testValue'] >=TT :
            rawCPlist[genename]= abs(valueinformation['ratio'])
    
    orderedrawCPlist = dict(sorted(rawCPlist.items(), key=lambda item: item[1],reverse=True))

    orderedgenelist={}
    index = 1
    for genename in orderedrawCPlist.keys():
        orderedgenelist[genename] = index
        index+=1


    return orderedgenelist

def findenrichmentscore(pathslist,orderedgenelist):
    scoredict = {}
    for pathsname,genelist in pathslist.items():
        scorelist=[]
        for gene in genelist:
            if gene in orderedgenelist.keys():
                scorelist.append(orderedgenelist[gene])
        if len(scorelist)>0:
            scoredict[pathsname] = mean(scorelist)
        else:
            scoredict[pathsname] = sys.maxsize

    
    return dict(sorted(scoredict.items(), key=lambda item: item[1]))



def writedicttofile(scoredict,outputfile):
    with open(outputfile,'w') as output:
        for pathname,pathdescore in scoredict.items():
            output.write(pathname+'\t'+str(pathdescore)+'\n')



def findtheenrichmentroutes(scorecollectiondict,topthresholds,sampleratio):

    hightestenrichmentroutelist={}
    outputenrichmentroutelist=[]
    for filename,routescorelist in scorecollectiondict.items():
        index = 0

        for routename,enrichscore in routescorelist.items():
            if index<topthresholds:
                if routename not in hightestenrichmentroutelist.keys():
                    hightestenrichmentroutelist[routename]=1
                else:
                    hightestenrichmentroutelist[routename]+=1
            index+=1

    for routename, topcontainsnumber in hightestenrichmentroutelist.items():
        if  (topcontainsnumber/len(scorecollectiondict))>=sampleratio:
            outputenrichmentroutelist.append(routename)
    
    return outputenrichmentroutelist










if __name__=="__main__":
    # load route path

    routepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_7-4/uniquepaths.txt"
    pathslist= loadroutes(routepath)
    outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/GSEA7-18/enrichmentscore_SKCM/"

    # ratiodir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/testHFD/"
    ratiodir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_Individual/"
    

    
    g = os.walk(ratiodir)  

    scorecollectionlist={}


    for path,dir_list,file_list in g:  
        for filename in file_list:
           ratiofile = os.path.join(path, filename) 
           CPlist = loadratiofile(ratiofile)
           orderedgenelist = orderDEgene(CPlist,TC=0,TT=0)
           scorelist = findenrichmentscore(pathslist,orderedgenelist)

           scorecollectionlist[filename]=scorelist
           outputfile = outputdir+filename
           writedicttofile(scorelist,outputfile)

    
    outputroutepath = findtheenrichmentroutes(scorecollectionlist,1500,0.8)

    with open(outputdir+'summary.txt','w') as summaryfile:
        for routename in  outputroutepath:
            summaryfile.write(routename+'\n')


