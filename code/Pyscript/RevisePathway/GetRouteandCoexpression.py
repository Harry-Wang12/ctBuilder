

import subprocess
import dataprocess
import CoexpressionDB
import sys
import os
import statistics 
from skrebate import ReliefF
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from multiprocessing import Process
from os import system



def getCoexpressionGene(genename,Coexpresslist):
    Coexpressionlist= []
    for genepair in Coexpresslist:
        if genepair["gene1"] == genename and not genepair["gene2"] in Coexpressionlist:
            Coexpressionlist.append(genepair["gene2"])
        elif genepair["gene2"] == genename and not genepair["gene1"] in Coexpressionlist:
            Coexpressionlist.append(genepair["gene1"])
    
    return Coexpressionlist
            





def generateRouteAndCoexpressionfileDown(datafileexpresionfilepath, PathwayRoutefilepath,CoexpressionDir,subgroupdir,subgroupnumber,accuvaluethreshold,importancethreshold):
    accuvaluethreshold = int(accuvaluethreshold*subgroupnumber)
    # l = []
    # hunman
    biogridfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/BiogridDB/BIOGRID-ORGANISM-Homo_sapiens-4.2.192.tab3.txt"
    Stringlinkfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Homo sampines/9606.protein.links.v11.0.txt"
    Stringinfofilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Homo sampines/9606.protein.info.v11.0.txt"
    # mouse
    # biogridfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/BiogridDB/BIOGRID-ORGANISM-Mus_musculus-4.2.192.tab3.txt"
    # Stringlinkfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Mus musculus/10090.protein.links.v11.0.txt"
    # Stringinfofilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Mus musculus/10090.protein.info.v11.0.txt"



    Coexpresslist1 = CoexpressionDB.loadfileBiogrid(biogridfilepath)
    Coexpresslist2 = CoexpressionDB.loadfileString(Stringlinkfilepath,Stringinfofilepath)
    Coexpresslist = CoexpressionDB.combinetwodb(Coexpresslist1,Coexpresslist2)


    # # # using php to generate Scoreing first
    # print('C:/wamp64/bin/php/php5.6.40/php.exe', 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/score_pathways_cli.php',datafileexpresionfilepath,PathwayRoutefilepath)

    # # sys.exit()
    # if you want output
    result = subprocess.run(
        ['C:/wamp64/bin/php/php5.6.40/php.exe', 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/score_pathways_cli.php',datafileexpresionfilepath,PathwayRoutefilepath],    # program and arguments
        stdout=subprocess.PIPE,  # capture stdout
        check=False               # raise exception if program fails
    )   

    # # sys.exit()

    expressiondata = dataprocess.loadfilematrix(datafileexpresionfilepath)


    with open(PathwayRoutefilepath,"r") as PathwayRoutefile:
        for line in PathwayRoutefile.readlines():
            linedata = line.strip().split("\t")
            if "~p2~" in linedata[1]:
                Secondline = linedata[1].replace("/","").split("~")
                Sourcelist = Secondline[3].split(",")
                Sourceid = linedata[0].split("~")[2]
                Targetlist = Secondline[5].split(",")
                Targetid = linedata[0].split("~")[3]
                # Contains= linedata[1].split("~p2~")[1].split(",")
                Coexpressionlist = []
                for genename in Sourcelist:
                    dbselectlist = getCoexpressionGene(genename,Coexpresslist)
                    Coexpressionlist=list(set(Coexpressionlist).union(set(dbselectlist)))
                                              
                with open(CoexpressionDir+linedata[1].split("~p2~")[0]+"_"+Secondline[3]+"_"+Sourceid+"_"+Targetid+".txt","w") as smallfile:
                    title = "Samplename"
                    for i in range(2,len(linedata)):
                        title+="\tS"+str(i-1)
                    smallfile.write(title+"\n")

                    p2line ="p2"
                    for i in range(2,len(linedata)):
                        p2line+="\t"+str(linedata[i]) 
                    smallfile.write(p2line+"\n")                
                    for genename in Coexpressionlist :
                        if genename in expressiondata.keys() and not genename in Sourcelist and not genename in Targetlist:
                            pline = genename
                            for value in expressiondata[genename]:
                                pline+="\t"+str(value) 
                            smallfile.write(pline+"\n")
                # do sub group
                dirname = linedata[1].split("~p2~")[0]+"_"+Secondline[3]+"_"+Sourceid+"_"+Targetid
                handleAroute(subgroupnumber,dirname,CoexpressionDir,subgroupdir,accuvaluethreshold,importancethreshold)
                


def handleAroute(subgroupnumber,dirname,CoexpressionDir,subgroupdir,accuvaluethreshold,importancethreshold):
    print(subgroupdir+dirname)
    dataprocess.randomdatasplit(subgroupnumber,CoexpressionDir+dirname+".txt",subgroupdir+dirname+"/")
    Reliffeatureselection(subgroupnumber,subgroupdir+dirname+"/","p2")    
    Findcorrectgene("P2",subgroupdir+dirname+"/result/",subgroupdir+dirname+"/analy/",accuvaluethreshold,isreg=False,ascend=True,toprate=importancethreshold)
    



def processRelif(index,datadir,targetgene):
    # print(index)      
    genetic_data = pd.read_csv(datadir+str(index)+'_test.txt',  sep='\t')
    
    features, labels = genetic_data.drop(targetgene, axis=1).values, genetic_data[targetgene].values
    # Make sure to compute the feature importance scores from only your training set
    X_train, X_test, y_train, y_test = train_test_split(features, labels)
    fs = ReliefF()
    fs.fit(X_train, y_train)
    with open(datadir+"/result/"+str(index)+"_result.txt","w") as resultfile:
        for feature_name, feature_score in zip(genetic_data.drop(targetgene, axis=1).columns,fs.feature_importances_):
            resultfile.write(feature_name+'\t'+targetgene+'\t'+str(feature_score)+"\n")


def Reliffeatureselection(totalsubgroup,datadir,targetgene):
    l=[]
    targetgene = targetgene.upper()
    for index in range(totalsubgroup):
        # processRelif(index,datadir,targetgene)
        p = Process(target = processRelif,args=(index,datadir,targetgene))
        p.start()
        l.append(p) 
    
    for p in l :
        p.join() 

        # print(index)      
        # genetic_data = pd.read_csv(datadir+str(index)+'_test.txt',  sep='\t')
        
        # features, labels = genetic_data.drop(targetgene, axis=1).values, genetic_data[targetgene].values
        # # Make sure to compute the feature importance scores from only your training set
        # X_train, X_test, y_train, y_test = train_test_split(features, labels)
        # fs = ReliefF()
        # fs.fit(X_train, y_train)
        # with open(datadir+"/result/"+str(index)+"_result.txt","w") as resultfile:
        #     for feature_name, feature_score in zip(genetic_data.drop(targetgene, axis=1).columns,fs.feature_importances_):
        #         resultfile.write(feature_name+'\t'+targetgene+'\t'+str(feature_score)+"\n")

def Findcorrectgene(findname,dirpath,resultpath,accuvaluethreshold,isreg=False,ascend=True,toprate=0.1):
    g = os.walk(dirpath)  
    resultslist=[]
    for path,dir_list,file_list in g:  
        for file_name in file_list:
            testresult = {}
            addlist = []
            with open(os.path.join(path, file_name),"r") as resultfile:
                for line in resultfile.readlines():
                    linedata = line.strip().split("\t")
                    TF = linedata[0].upper()

                    target = linedata[1].upper()
                    
                    weight = float(linedata[2].upper())
                    if isreg:
                        if TF == findname:
                            testresult[target] = weight
                    else:
                        if target == findname:
                            testresult[TF] = weight

            testresult={k: v for k, v in sorted(testresult.items(), key=lambda item: item[1],reverse=ascend)}
            accpectnum = toprate*len(testresult)
            for key,valueweight in testresult.items():
                accpectnum-=1
                if accpectnum>=0:
                    addlist.append(key)            
            resultslist.append(addlist)
    selectgene ={}
    for addedlist in resultslist:
        for genename in addedlist:
            if genename not in selectgene.keys():
                selectgene[genename]=1
            else:
                selectgene[genename]+=1
    if isreg:
        filename = findname+"_target.txt"
        wastefilename = findname+"_targetwaste.txt"
    else:
        filename = findname+"_regulator.txt"
        wastefilename = findname+"_regulatorwaste.txt"
    with open(resultpath+filename,"w")as resultouputfile:
        with open(resultpath+wastefilename,"w")as resultouputwastefile:
            for genenamekey,accuvalue in selectgene.items():
                if accuvalue>= accuvaluethreshold:
                    resultouputfile.write(genenamekey+"\t"+str(accuvalue)+"\n")
                    resultouputwastefile.write(genenamekey+"\t"+str(accuvalue)+"\n")
                else:
                    resultouputwastefile.write(genenamekey+"\t"+str(accuvalue)+"\n")



# if __name__ == '__main__':

    # importancethreshold = 0.1
    # subgouppercentage = 0.7
    # subgroupnumber = 10
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # outputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-27-2021/TCGA_SKCM_RouteScore.txt"
    # CoexpressionDir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/CoexpressionDir/"
    # subgroupdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/CoexpressionSubgroup/"
    # generateRouteAndCoexpressionfileDown(inputtable, outputtable,CoexpressionDir,subgroupdir,subgroupnumber,subgouppercentage,importancethreshold)
    
    # datadir= "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-18-2021/CoexpressionSubgroup/[8]Adherens junction/"
    # Reliffeatureselection(10,datadir,"P2")

    # os.system("shutdown -s -t  1")