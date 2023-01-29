
import os
import statistics 
import json



def loadfilematrix(filepath,transpose=False,title=False):
    returndict = {}

    if not transpose:
        with open(filepath,"r") as datafile:
            for line in datafile.readlines():
                if title:
                    title = False
                else:
                    linedate = line.strip().split("\t")
                    keyname = linedate[0]
                    del(linedate[0])
                    returndict[keyname] = linedate
    
    return returndict


def loadfilelog2dir(dirpath,genename=0,logratio=1,Cexpression=2,Texpression=3, CAVGexpression = False, CCZexpression = False,TAVGexpression = False, TCZexpression = False):
    # genename=0,log2=1,Cexpression=2,Texpression=3, CAVGexpression = 4, CCZexpression = 5,TAVGexpression = 6, TCZexpression = 7

    returnlist ={
        "totalgenlist":[],
        "detailcontent":[]
    }

    g = os.walk(dirpath)  

    for path,dir_list,file_list in g:  
        for file_name in file_list: 
            # print(len(file_name)) 
            filecontent = {
                "name": file_name,
                "log2": {},
                "genelist":[]
            }
            if not Cexpression==False:
                filecontent["Cexpression"]={}
            if not Texpression==False:
                filecontent["Texpression"]={}
            if not CAVGexpression==False:
                filecontent["CAVGexpression"]={}
            if not CCZexpression==False:
                filecontent["CCZexpression"]={}
            if not TAVGexpression==False:
                filecontent["TAVGexpression"]={}
            if not TCZexpression==False:
                filecontent["TCZexpression"]={}
                
            filepath = os.path.join(path, file_name) 
            header = True
            for line in open(filepath):
                if not line.startswith("!"):
                    if header:
                        header=False
                        continue
                    else:
                        linedata = line.strip().split("\t")                        
                        geneN = linedata[genename].upper()
                        
                        ratio = float(linedata[logratio])
                        if geneN not in filecontent["genelist"]:
                            filecontent["genelist"].append(geneN)

                        filecontent["log2"][geneN] = ratio                                               

                        if not Cexpression ==False:
                            if len(linedata)>Cexpression:
                                filecontent["Cexpression"][geneN]=float(linedata[Cexpression])
                        if not Texpression ==False:
                            if len(linedata)>Texpression:
                                filecontent["Texpression"][geneN]=float(linedata[Texpression])
                        if not CAVGexpression ==False:
                            if len(linedata)>CAVGexpression:
                                filecontent["CAVGexpression"][geneN]=float(linedata[CAVGexpression])
                        if not CCZexpression ==False:
                            if len(linedata)>CCZexpression:
                                filecontent["CCZexpression"][geneN]=float(linedata[CCZexpression])
                        if not TAVGexpression ==False:
                            if len(linedata)>TAVGexpression:
                                filecontent["TAVGexpression"][geneN]=float(linedata[TAVGexpression])
                        if not TCZexpression ==False:
                            if len(linedata)>TCZexpression:
                                filecontent["TCZexpression"][geneN]=float(linedata[TCZexpression])
            returnlist["detailcontent"].append(filecontent)
            if len(returnlist["totalgenlist"])==0:
                returnlist["totalgenlist"]=filecontent["genelist"]
            else:
                returnlist["totalgenlist"] = [i for i in returnlist["totalgenlist"] if i in filecontent["genelist"]]
    

    return returnlist

                    
def ConsistenceFitlerlog2UnDecided(regulatorkey,targetkey,datakeydict,threshold):

    # Undecide 
    # return relation or False if not meet the threshold 
    
    if regulatorkey in datakeydict.keys() and targetkey in datakeydict.keys():
        minlen = min([len(datakeydict[regulatorkey]),len(datakeydict[targetkey])])


        numActive = 0
        numInhibit = 0


        for i in range(minlen):
            if float(datakeydict[regulatorkey][i])*float(datakeydict[targetkey][i]) > 0:
                numActive+=1
            
            if float(datakeydict[regulatorkey][i])*float(datakeydict[targetkey][i]) < 0:
                numInhibit+=1
            
        if (numActive/minlen) >= threshold:
            return "activate"
        
        elif (numInhibit/minlen) >= threshold:
            return "inhibit"
        
        else:
            return False
    
    else:
        return False


def ConsistenceFitlerlog2Coexpress(genename,datadict,threshold):

    # Undecide 
    # Return True if it higher than threshold. Return False lower than threshold
    if genename in datadict["totalgenlist"]:
        numActive = 0
        numInhibit = 0
        for datacontent in datadict["detailcontent"]:
            if datacontent["log2"][genename] > 0:
                numActive+=1
            
            if datacontent["log2"][genename] < 0:
                numInhibit+=1

        if (numActive/len(datadict["detailcontent"])) >= threshold:
            return "activate"
        
        elif (numInhibit/len(datadict["detailcontent"])) >= threshold:
            return "inhibit"
        
        else:
            return False

    else:
        return False



    
def tranposmatrix(filepath,outputfilepath):
    originaldatadict = loadfilematrix(filepath)
    index=0
    with open(outputfilepath,"w") as outputfile:
        correctgenelist =[]
        samplenumber = 0
        for key,valuelist in originaldatadict.items():
            samplenumber=len(valuelist)
            correctgenelist.append(key)

        title="genename"
        for genename in correctgenelist:
            title+="\t"+genename
        
        outputfile.write(title+"\n")
        for i in range(samplenumber):
            index+=1
            pline = "Sample"+str(index)               
            for genename in correctgenelist:
                pline+="\t"+str(originaldatadict[genename][i])
            outputfile.write(pline+"\n")



    

def randomdatasplit(numberofgroup,datafilepath,outputdir,searchinggenelistfilepath=False,columnname=False,rowname = False,stdlarge=False):
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        os.makedirs(outputdir+"result/")
        os.makedirs(outputdir+"analy/")

    # searchinggenelist = []

    # with open(searchinggenelistfilepath,"r") as allgene:
    #     for line in allgene.readlines():
    #         linedata = line.strip().upper()
    #         searchinggenelist.append(linedata)

    resultdata=[]

    for i in range(numberofgroup):
        resultdata.append({})
  

    with open(datafilepath,"r") as datafile:
        title=True
        for line in datafile.readlines():
            if title:
                title = False
            else:
                linedata= line.strip().split("\t")
                
                genename = linedata[0].upper()
                # if genename in searchinggenelist:
                del(linedata[0])
                for i in range(numberofgroup):
                    resultdata[i][genename] = []
                for j in range(len(linedata)):
                    resultdata[j%numberofgroup][genename].append(float(linedata[j]))

    # write
    index =0
    for record in resultdata:        
        genelist=[]
        with open(outputdir+str(index)+"_test.txt","w") as recorddata:
            correctgenelist =[]
            samplenumber = 0
            for key,valuelist in record.items():
                samplenumber=len(valuelist)
                if stdlarge:
                    if statistics.stdev(valuelist):
                        correctgenelist.append(key)
                else:
                    correctgenelist.append(key)

            title=""
            for genename in correctgenelist:
                title+=genename+"\t"
            
            recorddata.write(title[:-1]+"\n")

            for i in range(samplenumber):
                pline = ""               
                for genename in correctgenelist:
                    pline+=str(record[genename][i])+"\t"
                recorddata.write(pline[:-1]+"\n")
        index+=1


   
def getpathwaygenes(pathwaydir,outputdir):

    g = os.walk(pathwaydir)

    for path,dir_list,file_list in g:  
        for file_name in file_list: 
            with open(os.path.join(path, file_name),"r") as f:
                pathway = json.load(f)
                with open(outputdir+"Allgene_"+file_name.split(".")[0]+".txt","w") as allgenefile:
                    for node in pathway["elements"]["nodes"]:
                        if len(node["data"]["name"])>0 and not node["data"]["name"].startswith("n") and not node["data"]["name"].istitle():
                            allgenefile.write(node["data"]["name"]+"\n")
                       

    
   



           
if __name__=="__main__":

    pathwaydir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/pathways/"
    outputdir   = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/pathways/Keggpathwaygene/"
    getpathwaygenes(pathwaydir,outputdir)



  