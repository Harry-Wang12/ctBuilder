import statistics
from scipy import stats
import numpy as np
import os

def loadrandomfile(randomfilepath):
    index = 0 
    returndict = {}
    with open(randomfilepath,"r") as randomfile:
        for line in randomfile.readlines():
            index+=1
            if index%1000 ==0:
                print("load random dict "+str(index))
            linedata=line.strip().split("\t")
            target = linedata[0]
            source = linedata[1]
            del(linedata[0])
            del(linedata[0])
            if not target  in returndict.keys():
                returndict[target]={
                   source:[] 
                }
            elif not source  in returndict[target].keys():
                returndict[target][source]=[]
            # print(target,source)
            for weight in linedata:
                returndict[target][source].append(float(weight))
    
    return returndict

def reviseresultfile(resultfilepath,normaldict):    
    returnlines=[]
    with open(resultfilepath,"r") as resultfile:
        for line in resultfile.readlines():
            linedata=line.strip().split("\t")
            target = linedata[0]
            source = linedata[1]           
            weight = float(linedata[2])
            print(target,source)
            if target in normaldict.keys():
                if source in normaldict[target].keys():
                    ttestresult = stats.ttest_1samp(normaldict[target][source],weight)
                else:
                    ttestresult = stats.ttest_1samp( [0]*100,weight)
            else:
                    ttestresult = stats.ttest_1samp([0]*100,weight)

            
            pvalue = ttestresult.pvalue 

            pline = target+'\t'+source+'\t'+str(weight)+'\t'+str(pvalue)+'\n'

            returnlines.append(pline)
    
    with open(resultfilepath,"w") as resultfile:
        for line in returnlines:
            resultfile.write(line)



def removeunwantedgefromfile(resultfilepath,importancethreshold=1,pvaluethreshold = 0.05):

    resultlist=[]

    with open(resultfilepath,'r') as resultfile:
        
        for line in resultfile.readlines():
            
            linedata = line.strip().split("\t")
            if len(linedata)>3 and not linedata[3]=='nan':
                edgeinfo = {
                    'target':linedata[0],
                    'source':linedata[1],
                    'importance':float(linedata[2]),
                    'pvalue':float(linedata[3])
                }
            else:
                edgeinfo = {
                    'target':linedata[0],
                    'source':linedata[1],
                    'importance':float(linedata[2]),
                    'pvalue':'nan'
                }
            
            resultlist.append(edgeinfo)

    return removeunwantedgefromlist(resultlist,importancethreshold,pvaluethreshold)


def removeunwantedgefromlist(edgelist,importancethreshold=1,pvaluethreshold = 0.05):
    resultlist=[]
    
    for edgeinfo in edgelist:
        if edgeinfo['pvalue'] == 'nan' or edgeinfo['pvalue'] < pvaluethreshold:
            if edgeinfo['importance'] >=importancethreshold:
                resultlist.append(edgeinfo)
    
    return resultlist



def writeedgetofile(outputfilepath,edgelist):
    with open(outputfilepath,'w') as outputfile:
        for edgeinfo in edgelist:
            outputfile.write(edgeinfo['target']+'\t'+edgeinfo['source']+'\t'+str(edgeinfo['importance'])+'\t'+str(edgeinfo['pvalue'])+'\n')



        
     
def findoutCommonedge(dirpath,outputfilepath,remainrate=0.9,importancerate=1.4):
    
    fileremain={}
    filenumber = 0
    
    files= os.listdir(dirpath)
    for file in files: 
        if not os.path.isdir(file): 
            refinefilepath = dirpath+file  
            filenumber+=1          
            # the first file
            with open(refinefilepath,'r') as refinefile:                
                for line in refinefile.readlines():
                    linedata = line.strip().split('\t')
                    target = linedata[0]
                    source = linedata[1]
                    importance = float(linedata[2])
                    if importance>=importancerate:
                        if not target in fileremain.keys():
                            fileremain[target]={}
                        if not source in fileremain[target].keys():
                            fileremain[target][source]=1
                        else:
                            fileremain[target][source]+=1

    with open(outputfilepath,'w') as outputfile:
        for target,sources in fileremain.items():            
            for source,aggreenumber in sources.items():
                if aggreenumber/filenumber>=remainrate:
                    outputfile.write(source+'\t'+target+'\n')


     
def findoutCommonedgerate(dirpath,outputfilepath,remainrate=0.9,importancerate=0.25,minscore =1):
    
    fileremain={}
    filenumber = 0
    files= os.listdir(dirpath)
    for file in files: 
        if not os.path.isdir(file): 
            refinefilepath = dirpath+file  
            filenumber+=1          
            # the first file
            filecontentdict = {}
            with open(refinefilepath,'r') as refinefile:                
                for line in refinefile.readlines():
                    linedata = line.strip().split('\t')
                    target = linedata[0]
                    source = linedata[1]
                    importance = float(linedata[2])                  

                    if not target in filecontentdict.keys():
                        filecontentdict[target]={}
                    if not source in filecontentdict[target].keys() and importance>minscore:
                        filecontentdict[target][source]=importance

            for target , sourcedict in filecontentdict.items():
                    orderedsourcelist = dict(sorted(sourcedict.items(), key=lambda item: item[1],reverse=True)) 
                    selectedvalue = int(len(orderedsourcelist)*importancerate)+1
                    index =0
                    for sourcename,importancevalue in orderedsourcelist.items():
                        if index<=selectedvalue:
                            if not target in fileremain.keys():
                                fileremain[target]={}
                            if not sourcename in fileremain[target].keys():
                                fileremain[target][sourcename]=1
                                index+=1
                            else:
                                fileremain[target][sourcename]+=1
                                index+1
                            

                  

    with open(outputfilepath,'w') as outputfile:
        for target,sources in fileremain.items():            
            for source,aggreenumber in sources.items():
                if aggreenumber/filenumber>=remainrate :
                    outputfile.write(source+'\t'+target+'\n')

    




def mappwithalledge(filenamepathlist,outputfilepath):

    resultdict={}
    originalfilenamelist = []
    for filenamepath in filenamepathlist :
        with open(filenamepath,'r')as filename:
            originalfilename  = filenamepath.split('/')[-1]
            originalfilenamelist.append(originalfilename)
            for line in filename.readlines():
                linedata= line.strip().split('\t')
                target = linedata[0]
                source = linedata[1]
                if not target in resultdict.keys():
                    resultdict[target]={
                        source:[originalfilename]
                    }
                if not source in resultdict[target].keys():
                    resultdict[target][source]=[originalfilename]
                if not originalfilename in resultdict[target][source]:
                    resultdict[target][source].append(originalfilename)

    with open(outputfilepath,'w') as outputfile:
        title = "target\tsource"
        for filename in originalfilenamelist:
            title+="\t"+filename        
        outputfile.write(title+"\n")

        for target,sourcelist in resultdict.items():
            for source, filenames in sourcelist.items():
                
                pline= target+'\t'+source
                for filename in originalfilenamelist:
                    if filename in filenames:
                        pline+='\t1'
                    else:
                        pline+='\t0'
                outputfile.write(pline+"\n")


    



def findcommonedge(filepath1, filepath2,outputfilepath):
    edgelist1={}
    with open(filepath1,'r') as refinefile:                
                for line in refinefile.readlines():
                    linedata = line.strip().split('\t')
                    target = linedata[0]
                    source = linedata[1]
                               

                    if not target in edgelist1.keys():
                        edgelist1[target]=[]
                    if not source in edgelist1[target]:
                        edgelist1[target].append(source)
    

    commonedge = {}
    with open(filepath2,'r') as refinefile:                
                for line in refinefile.readlines():
                    linedata = line.strip().split('\t')
                    target = linedata[0]
                    source = linedata[1]                               

                    
                    if target in edgelist1.keys():
                       if  source in edgelist1[target]:
                            if not target in commonedge.keys():
                                commonedge[target]=[]
                            if not source in commonedge[target]:
                                commonedge[target].append(source)

    with open(outputfilepath,'w') as outputfile:
        for target,sources in commonedge.items():            
            for source in sources:
                outputfile.write(source+'\t'+target+'\n')

    return commonedge
                    
    




            
           
if __name__=="__main__":


    outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/remainlink/"

    # controldir = outputdir+"kmeans_control/"
    testdir =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-17/singlecellAD/"

    remainrate=0.8
    importancerate=0.1
    # controloutputfilepath = outputdir+"controlremain.txt"
    testoutputfilepath = outputdir+"testremain_AD_"+str(remainrate)+"_"+str(importancerate)+".txt"
    # randomputfilepath = outputdir+"testremain_random_0.8_0.1.txt"

    # outputfilepath = outputdir+"common_test_random_0.8_0.1.txt"
    # findoutCommonedge(controldir,controloutputfilepath,remainrate=0.55,importancerate=1.25)
    # findoutCommonedge(testdir,testoutputfilepath,remainrate=0.55,importancerate=1.25)


    testoutputfilepath = outputdir+"testremain_random_0.8_0.1.txt"
    randomputfilepath = outputdir+"testremain_AD_0.8_0.1.txt"
    outputfilepath =  outputdir+"common_random_AD.txt"
    findcommonedge(testoutputfilepath, randomputfilepath,outputfilepath)
            
    # findoutCommonedgerate(controldir,controloutputfilepath,remainrate=0.9,importancerate=0.05)
    # findoutCommonedgerate(testdir,testoutputfilepath,remainrate=remainrate,importancerate=importancerate)


            
            


