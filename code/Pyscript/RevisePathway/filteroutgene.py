import torch
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pyplot as plt
import numpy as np
import dataprocess
import os
from GENIE3 import *
import pandas as pd
import numpy as np
from sklearn.pipeline import make_pipeline
from skrebate import ReliefF
from sklearn.model_selection import train_test_split
import statistics 


def normalization(data):
    _range = np.max(data) - np.min(data)
    return (data - np.min(data)) / _range


class Net(nn.Module):
    # 6*3
    def __init__(self,n_input,n_hidden,n_output):
        super(Net,self).__init__()
        # self.hidden1 = nn.Linear(n_input,n_hidden)
        self.predict = nn.Linear(n_input,n_output)
    def forward(self,input):
        # out = self.hidden1(input)
        # out =F.sigmoid(out)
        out = self.predict(input)
        out =F.sigmoid(out)

        return out


def Method1(upcohortdirTrain,downcohortdirTrain,upcohortdirTEST,downcohortdirTEST,allgenelistfile,traininggenelist,remainedgenelistfile,lr,epchos,outputfile,TFgene="IRF1"):
    upcohorttrain = dataprocess.loadfilelog2dir(upcohortdirTrain)
    downcohorttrain = dataprocess.loadfilelog2dir(downcohortdirTrain)

    upcohorttest = dataprocess.loadfilelog2dir(upcohortdirTEST)
    downcohorttest = dataprocess.loadfilelog2dir(downcohortdirTEST)

    allgenelist = []
    correctgenelist=[]
    modifiletrain=[]
    with open(allgenelistfile,"r") as allgene:
        for line in allgene.readlines():
            linedata = line.strip().upper()
            allgenelist.append(linedata)
    
    with open(remainedgenelistfile,"r") as allgene:
        for line in allgene.readlines():
            linedata = line.strip().upper()
            correctgenelist.append(linedata)
    
    with open(traininggenelist,"r") as allgene:
        for line in allgene.readlines():
            linedata = line.strip().upper()
            modifiletrain.append(linedata)
    
    
    trainingX= []
    trainingY = []
    testX= []
    testY = []

    minlentrain = min([len(upcohorttrain["detailcontent"]),len(downcohorttrain["detailcontent"])])
    minlentest = min([len(upcohorttest["detailcontent"]),len(downcohorttest["detailcontent"])])

    for genename in allgenelist:
        if genename in upcohorttrain["totalgenlist"] and genename in downcohorttrain["totalgenlist"] and genename in modifiletrain:
            for i in range(minlentrain):
                trainingrecord = [1,upcohorttrain["detailcontent"][i]["log2"][genename],upcohorttrain["detailcontent"][i]["log2"][TFgene],0,downcohorttrain["detailcontent"][i]["log2"][genename],downcohorttrain["detailcontent"][i]["log2"][TFgene]]
                trainingX.append(trainingrecord)
                
                if genename in correctgenelist:
                    trainingY.append([1])
                   
                else:
                    trainingY.append([0])
        if genename in upcohorttest["totalgenlist"] and genename in downcohorttest["totalgenlist"]:            
            for j in range(minlentest):
                testgrecord = [1,upcohorttest["detailcontent"][j]["log2"][genename],upcohorttest["detailcontent"][j]["log2"][TFgene],0,downcohorttest["detailcontent"][j]["log2"][genename],downcohorttest["detailcontent"][j]["log2"][TFgene]]
                testX.append(testgrecord)
                if genename in correctgenelist:
                    testY.append([1])
                else:
                    testY.append([0])
   
    net = Net(6,3,1)
    optimizer = torch.optim.SGD(net.parameters(),lr = lr)
    loss_func = torch.nn.MSELoss()
    trainingX = torch.from_numpy(np.array(trainingX)).float()
    trainingY = torch.from_numpy(np.array(trainingY)).float()

    traintotal=[]
    testtotal=[]
    for t in range(epchos):
        trainingloss=[]
        for i in range(len(trainingX)):
            net.train()
            X = torch.from_numpy(np.array(trainingX[i])).float()
            Y = torch.from_numpy(np.array(trainingY[i])).float()
            prediction = net(X)
            loss = loss_func(prediction,Y)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            trainingloss.append(loss)
        
        traintotal.append(sum(trainingloss) / len(trainingloss))
        with torch.no_grad(): 
            testloss=[]
            net.eval() 
            for j in range(len(testX)): 
                tX = torch.from_numpy(np.array(testX[j])).float()
                tY = torch.from_numpy(np.array(testY[j])).float() 
                testloss=[]
                predictionT = net(tX)
                tloss = loss_func(predictionT ,tY)
                # print(predictionT)
                testloss.append(tloss)
            print("epoch:",t,"training loss: ", sum(trainingloss) / len(trainingloss),"test loss: ", sum(testloss) / len(testloss))
            testtotal.append(sum(testloss) / len(testloss))
    
    with open(outputfile,"w") as resultcomparefile:
        resultcomparefile.write("Genename\tUp\ttargetlog2\tTFlog2\tDown\ttargetlog2\tTFlog2\tPredict\tLabel\n")
        for j in range(minlentest):
            for genename in allgenelist:           
                if genename in upcohorttest["totalgenlist"] and genename in downcohorttest["totalgenlist"]: 
                    testgrecord = [1,upcohorttest["detailcontent"][j]["log2"][genename],upcohorttest["detailcontent"][j]["log2"][TFgene],0,downcohorttest["detailcontent"][j]["log2"][genename],downcohorttest["detailcontent"][j]["log2"][TFgene]]
                    tX = torch.from_numpy(np.array(testgrecord)).float()
                    predictionT = net(tX)
                    if predictionT >= 0.5:
                        prediction = "1"
                    else:
                        prediction = "0"

                    if genename in correctgenelist:
                        label = "1"
                    else:
                        label = "0"

                    resultcomparefile.write(genename+"\tUp\t"+str(upcohorttest["detailcontent"][j]["log2"][genename])+"\t"+str(upcohorttest["detailcontent"][j]["log2"][TFgene])+"\tDown\t"+str(downcohorttest["detailcontent"][j]["log2"][genename])+"\t"+str(downcohorttest["detailcontent"][j]["log2"][TFgene])+"\t"+prediction+"\t"+label+"\n")
                
    X=range(epchos)
    
    l1=plt.plot(X,traintotal,'r--',label='train')
    l2=plt.plot(X,testtotal,'g--',label='test')
    
    plt.plot(X,traintotal,'ro-',X,testtotal,'g+-')

    plt.xlabel('iter')
    plt.ylabel('loss')
    plt.legend()
    plt.show()






# network related
def randomdatasplit(numberofgroup,datafilepath,searchinggenelistfilepath,outputdir,columnname=False,rowname = False,stdlarge=False):
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    searchinggenelist = []

    with open(searchinggenelistfilepath,"r") as allgene:
        for line in allgene.readlines():
            linedata = line.strip().upper()
            searchinggenelist.append(linedata)

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
                if genename in searchinggenelist:
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
                        
                        # pline = ""
                        # for value in valuelist:
                        #     pline+=str(value)+"\t"
                        # pline=pline[:-1]+"\n"
                        # recorddata.write(pline)
                        # genelist.append(key)
        
        # with open(outputdir+str(index)+"_genename.txt","w") as recorddata:
        #     pline = ""
        #     for key in genelist:
        #         pline+=key+"\t"
        #         pline=pline[:-1]+"\n"
        #     recorddata.write(pline)
        index+=1


   

    # print(resultdata)

   

def usingGENIN3getresult(dir,groupnumber):

    if not os.path.exists(dir+"result/"):
        os.makedirs(dir+"result/")

    for i in range(groupnumber):
        
        gene_names=[]
        with open(dir+str(i)+"_genename.txt","r") as f:
            for line in f.readlines():
                gene_names.append(line.strip()) 
            
        data = np.loadtxt(dir+str(i)+"_test.txt")
        data = data.transpose()
        # data = normalization(data)
        VIM = GENIE3(data,nthreads=8)
        get_link_list(VIM,gene_names=gene_names,file_name=dir+"result/"+str(i)+"_result.txt")
        





def Findcorrectgene(findname,dirpath,resultpath,isreg,ascend=True,toprate=0.1):
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
    else:
        filename = findname+"_regulator.txt"
    with open(resultpath+filename,"w")as resultouputfile:
        for genenamekey,accuvalue in selectgene.items():
            resultouputfile.write(genenamekey+"\t"+str(accuvalue)+"\n")





def Reliffeatureselection():
    # for index in range(10):
    index = 1
    genetic_data = pd.read_csv('C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/randomselect/test7_ccl2/'+str(index)+'_test.txt',  sep='\t')

    featurelist = ['NME2','TBX1','KIN','FOXO']
    refinefeaturelist = []

    for gene in featurelist:
        if gene in genetic_data.index:
            refinefeaturelist.append(gene)



    features1, labels = genetic_data.drop('CCL2', axis=1).values, genetic_data['CCL2'].values

    features2, labels = genetic_data[refinefeaturelist].values, genetic_data['CCL2'].values

    # Make sure to compute the feature importance scores from only your training set
    X_train, X_test, y_train, y_test = train_test_split(features2, labels)

    fs = ReliefF()
    fs.fit(X_train, y_train)

    with open("C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/randomselect/test7_ccl2/result/"+str(index)+"_result.txt","w") as resultfile:
        for feature_name, feature_score in zip(featurelist,fs.feature_importances_):
            resultfile.write(feature_name+'\tCCL2\t'+str(feature_score)+"\n")
            
            


    # print(resultdata)







Reliffeatureselection()

# Findcorrectgene(findname,dirpath,resultpath,isreg,toprate=0.1)





        
   












    







