import TFtargetnetworkrelated
import os
from itertools import permutations
import Graphrelated
import subprocess
import sys
import os
from multiprocessing import Process
import itertools

def getgenefromroute(outputdir,routefilepath):

    routegenelistdict={}

    with open(routefilepath,'r') as routefile:
        isheader = True
        for line in routefile.readlines():
            if isheader:
                isheader=False
            else:
                linedata = line.strip().replace('"','').split('\t')
                genelist = linedata[1].split("~")[5].split(",")
                routegenelistdict[linedata[0]]=genelist
    
    for routename,genelist in routegenelistdict.items():
        outputfilepath = outputdir+routename+'.txt'
        with open(outputfilepath,'w') as outputfile:
            for gene in genelist:
                outputfile.write(gene+"\n")


def loadtotalfile(allsamplefilepath):
    allsamples=[]

    with open(allsamplefilepath,'r') as allsamplefile:
        istitle=True

        for line in allsamplefile.readlines():
            linedata = line.strip().split("\t")
            if istitle:
                del(linedata[0])
                for sample in linedata:
                    allsamples.append({})
                istitle=False
            else:
                
                genename = linedata[0].upper()
                # print(genename)
                del(linedata[0])
                for i in range(len(linedata)):
                    allsamples[i][genename] = float(linedata[i])
    
    return allsamples


# def getweight(routepath):
#     n = len(routepath)-1
#     d =  


# def isconsistnce(gene1log,gene2log,gene1Z,gene2Z,pattern):
    

def findallgenecombine(routelist1,routelist2):

    combinelist = []
    for upgene in routelist1:
        for downgene in routelist2:
            if not upgene in routelist2 and not downgene in routelist1:
                combinelist.append([upgene,downgene])

    return combinelist

def readgenelist(filepath):
    genelist = []
    with open(filepath,'r') as genefile:
        for line in genefile.readlines():
            genelist.append(line.strip())
    
    return genelist

 
def getallgenelistfilename(dirpath):

    filelist=[]
    g = os.walk(dirpath)  

    for path,dir_list,file_list in g:  
        for file_name in file_list:  
            filelist.append(os.path.join(path, file_name))
    
    return filelist


def writepathtofile(pathlists,outputfilepath):

    with open(outputfilepath,'w') as  outputfile:
        for path in pathlists:
            splitstring='\t'
            outputfile.write(splitstring.join(path)+'\n')
                 

def findpathandwrite(graphdict,start,end,pathlength,outputfilepath,geneset1,geneset2):
    
    
    print("Start find route from "+start+" to "+end)       
    allpaths = Graphrelated.Check(graphdict,start,end,pathlength,geneset1,geneset2)
    
    writepathtofile(allpaths,outputfilepath)
    if len(allpaths)>0:
        print("Finish route from "+start+" to "+end +". "+str(len(allpaths))+" paths are found")
    return allpaths


def scoring(pathlist,patterns,logsample,Zsample=False):
	
    value = 0
    for i in range(len(patterns)):
        gene1 = pathlist[i]
        gene2 = pathlist[i+1]
        pattern=patterns[i]
        # Zscorestart = Zsample[gene1]
        # Zscoreend = Zsample[gene2]
        if gene1 in logsample.keys():
            logstart = logsample[gene1]
        else:
            logstart = 0

        if gene2 in logsample.keys():
            logend = logsample[gene2]
        else:
            logend = 0	

        if isconsistent(pattern,logstart,logend):
            value+=1

    return value/len(patterns)
	
def getupperstreampattern(edgelist,path):
    # print(path)
    # if path[0]=='CD14' and path[-2]=='MAPK14' and path[-1]=='FOS'  and path[-3]=='MAP2K3'and path[-4]=='MAP3K7' and path[1]=='TLR4' and path[2]=='TICAM2' and len(path)==10:
    #     print("11")
    patternlists=[]
    finalpatterns = []
    num2 =0
   
    
    for i in range(len(path)-1):
        if edgelist[path[i]][path[i+1]]['edge_type'] == 'activate':
            patternlists.append(1)
        elif edgelist[path[i]][path[i+1]]['edge_type'] == 'inhibit':
            patternlists.append(-1)
        elif edgelist[path[i]][path[i+1]]['edge_type'] == 'line':
            patternlists.append(0)
            num2+=1

    if num2> 0:                 
        combinations = itertools.product([-1,1], repeat=num2)
        
        for combination in combinations:
            pattern = []
            index = 0
            for i in range(len(path)-1):
                if patternlists[i]==1:
                    pattern.append(1)
                elif patternlists[i]==-1:
                    pattern.append(-1)
                elif patternlists[i]==0:
                    pattern.append(combination[index])
                    index +=1

            finalpatterns.append(pattern)
    else:
        pattern = []
        for i in range(len(path)-1):
            if patternlists[i]==1:
                pattern.append(1)
            elif patternlists[i]==-1:
                pattern.append(-1)
            # elif patternlists[i]==0:
            #     pattern.append(combination[index])
            #     index +=1

        finalpatterns.append(pattern)

    
    return finalpatterns

def getrealtion(relationfilepath):

    patternlists={}
    with open(relationfilepath,"r") as graphfile:
        for line in graphfile.readlines():
            linedata= line.strip().split("\t")
            tfgenename = linedata[0]
            targetgenename = linedata[1]
            relation = linedata[2]
            if not tfgenename in patternlists.keys():
                patternlists[tfgenename]={}

            if not targetgenename in patternlists[tfgenename].keys():
                patternlists[tfgenename][targetgenename]=relation

    return patternlists

def getmidperstreampattern(filepatternlists,path):

    patternlists=[]
    
    finalpatterns = []
    num2 =0

    for i in range(len(path)-1):
        if filepatternlists[path[i]][path[i+1]] == 'Activation':
            patternlists.append(1)
        if filepatternlists[path[i]][path[i+1]] == 'Repression':
            patternlists.append(-1)
        if filepatternlists[path[i]][path[i+1]] == 'Unknown':
            patternlists.append(0)
            num2+=1
    
    combinations = itertools.product([-1,1], repeat=num2)

    for combination in combinations:
        pattern = []
        index = 0
        for i in range(len(path)-1):
            if patternlists[i]==1:
                pattern.append(1)
            if patternlists[i]==-1:
                pattern.append(-1)
            if patternlists[i]==0:
                pattern.append(combination[index])
                index +=1

        finalpatterns.append(pattern)
    
    return finalpatterns

            


def isconsistent(pattern,logstart,logend, Zscorestart=0,Zscoreend=0):
	
	if pattern >0:
		if logstart * logend >0:
			return True
		elif Zscorestart * Zscoreend >0:
			return True
		else:
			return False
	
	if pattern <0:
		if logstart * logend <0:
			return True
		elif Zscorestart * Zscoreend <0:
			return True
		else:
			return False
	
	if pattern ==0:
		if logstart * logend >0:
			return True
		elif Zscorestart * Zscoreend >0:
			return True
		else:
			return False
			

# def getallgenefrompathwayfile()

# def 


def getpathfromfile(filepath):

    paths =[]
    with open(filepath,'r') as filex:
        for line in filex.readlines():
            paths.append(line.strip().split('\t'))

    return paths

# if __name__=="__main__":
    








