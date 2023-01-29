
import random
import os
import copy

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

def loadtotalfilewithname(allsamplefilepath):
    
    allsamples={}
    namelist = []
    with open(allsamplefilepath,'r') as allsamplefile:
        istitle=True
        for line in allsamplefile.readlines():
            linedata = line.strip().split("\t")
            if istitle:
                del(linedata[0])
                namelist = linedata
                for sample in linedata:
                    allsamples[sample]={}
                istitle=False
            else:
                
                genename = linedata[0].upper()
                # print(genename)
                del(linedata[0])
                for i in range(len(linedata)):
                    allsamples[namelist[i]][genename] = float(linedata[i])
    
    return allsamples


def seperateall(numberofgroup,allsamples,israndom):

    
    if israndom:
        random.shuffle(allsamples.keys())

    seperatedsamples=[]

    for i in range(numberofgroup):
        seperatedsamples.append({})

    index =-1
    for sample in allsamples:
        if index==numberofgroup-1:
            index =0
        else:
            index+=1
            
        seperatedsamples[index].append(sample)

    return seperatedsamples



def seperateallwithname(numberofgroup,allsamples,israndom):

    randomnamelist = []
    for samplename in allsamples.keys():
        randomnamelist.append(samplename)
    if israndom:
       
        random.shuffle(randomnamelist)
        

    seperatedsamples=[]

    for i in range(numberofgroup):
        seperatedsamples.append({})

    index =-1
    for sample in randomnamelist:

        if index==numberofgroup-1:
            index =0
        else:
            index+=1
            
        seperatedsamples[index][sample]= allsamples[sample]      

    return seperatedsamples


def splitmatrixwithtype(allsamplefilepath,typefilepath,samplecol,typecol,outputdir,needtype=False):
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    allsamples={}
    titleindex=[]
    genelist =[]

    with open(allsamplefilepath,'r') as allsamplefile:
        istitle=True

        for line in allsamplefile.readlines():
            linedata = line.strip().split("\t")
            if istitle:
                del(linedata[0])
                for sample in linedata:
                    allsamples[sample]={}
                    titleindex.append(sample)
                istitle=False
            else:                
                genename = linedata[0].upper()
                if not genename in genelist:
                    genelist.append(genename) 
                # print(genename)
                del(linedata[0])
                for i in range(len(linedata)):
                    allsamples[titleindex[i]][genename] = float(linedata[i])
    

    typedict ={}

    with open(typefilepath,'r') as typefile:
        istitle=True

        for line in typefile.readlines():
            linedata = line.strip().split("\t")
            
            if istitle:
                istitle=False
            else:   
                if len(linedata)>  typecol and len(linedata)>  samplecol :    
                    sampletype = linedata[typecol-1]
                    samplename = linedata[samplecol-1]
                    if sampletype not in typedict.keys():
                        typedict[sampletype] = [samplename]
                    else:
                        typedict[sampletype].append(samplename) 
        

    if needtype:
        outputfilepath = outputdir+needtype+".txt"
    # specific type
        allsamplenames = typedict[needtype]
        with open(outputfilepath,"w") as outputfile:
            titleline = "Gene"
            for samplename in allsamplenames :
                titleline+="\t"+samplename
            outputfile.write(titleline+'\n')

            for genename in genelist:
                pline = genename
                for samplename in allsamplenames:
                    pline+="\t"+str(allsamples[samplename][genename])
                outputfile.write(pline+'\n')
    
    else:
        for needtype in typedict.keys():
            outputfilepath = outputdir+needtype+".txt"
    # specific type
            allsamplenames = typedict[needtype]
            with open(outputfilepath,"w") as outputfile:
                titleline = "Gene"
                for samplename in allsamplenames :
                    titleline+="\t"+samplename
                outputfile.write(titleline+'\n')

                for genename in genelist:
                    pline = genename
                    print(genename)
                    for samplename in allsamplenames:
                        if samplename in allsamples.keys():
                            # print(samplename)
                            # print(genename)
                            pline+="\t"+str(allsamples[samplename][genename])
                    outputfile.write(pline+'\n')

def randomselectedsamples(selectednumber,allsamples):
    return random.sample(allsamples,selectednumber)


def splistmatrixandlabelwithname(numberofgroup,datasamples,labelsamples,israndom):




    randomnamelist = []
    for samplename in datasamples.keys():
        randomnamelist.append(samplename)
    if israndom:
       
        random.shuffle(randomnamelist)
        

    seperateddatasamples=[]

    for i in range(numberofgroup):
        seperateddatasamples.append({})

    index =-1
    for sample in randomnamelist:

        if index==numberofgroup-1:
            index =0
        else:
            index+=1
            
        seperateddatasamples[index][sample]= datasamples[sample]   


    seperatedlabelsamples=[]

    for i in range(numberofgroup):
        seperatedlabelsamples.append({})

    index =-1
    for sample in randomnamelist:

        if index==numberofgroup-1:
            index =0
        else:
            index+=1
            
        seperatedlabelsamples[index][sample]= labelsamples[sample]      

    return seperateddatasamples,seperatedlabelsamples
    



def writesampletofilename(allsamples,outputfile):

    samplenames=[]
    for samplename in allsamples.keys():
        samplenames.append(samplename)

    
    # # getallgene
    allgenes = allsamples[samplenames[0]].keys()

    # writting
    with open(outputfile,'w') as output:
        title = "genename"
        for samplename in samplenames:
            title+="\t"+samplename
        title = title+"\n"
        output.write(title)

        for genename in allgenes:
            pline = genename
            for samplename, samplegenevalues in allsamples.items():
                
                pline+="\t"+str(samplegenevalues[genename])

            pline = pline+"\n"
            output.write(pline)



def writesampletofile(allsamples,outputfile):
    # getallgene
    allgene = allsamples[0].keys()
    # writting
    with open(outputfile,'w') as output:
        title = ""
        for genename in allgene:
            title+=genename.upper()+"\t"
        title = title[:-1]+"\n"
        output.write(title)

        for sample in allsamples:
            pline = ""
            for genename in allgene:
                pline+=sample[genename.upper()]+"\t"
            pline = pline[:-1]+"\n"
            output.write(pline)



if __name__=="__main__":



    # kmeanscontrolsamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE169445_control_DC_combined.txt"
    # kmeanscontrollabelsamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeans_GSE169445_control_DC_combined.txt"
    # kmeanstestdatasamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE169445_Nash_DC_combined.txt"
    # kmeanstestlabelsamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeans_GSE169445_Nash_DC_combined.txt"

    kmeansrandomdatasamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/combineAD.txt"
    kmeansrandomlabelsamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/combineAD.txt"

    # controlallsamples = loadtotalfilewithname(kmeanscontrolsamples)

    # controlallsplitsamples = seperateallwithname(5,controlallsamples,True)


    testdatasamples = loadtotalfilewithname(kmeansrandomdatasamples)
    testlabelsamples = loadtotalfilewithname(kmeansrandomlabelsamples)

    testdatasplitsamples,testlabelsplitsamples = splistmatrixandlabelwithname(5,testdatasamples,testlabelsamples,True)

    index =1

    for testallsample in testdatasplitsamples:
        outputfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_"+str(index)+".txt"
        writesampletofilename(testallsample,outputfile)
        index+=1


    index =1

    for testallsample in testlabelsplitsamples:
        outputfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_label"+str(index)+".txt"
        writesampletofilename(testallsample,outputfile)
        index+=1


    # testallsplitsamples = seperateallwithname(5,testallsamples,True)

    # index =1

    # for controlallsample in controlallsplitsamples:
    #     outputfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed/kmeans_GSE169445_control_DC_combined"+str(index)+".txt"

    #     writesampletofilename(controlallsample,outputfile)
    #     index+=1
    

    # index =1

    # for testallsample in testallsplitsamples:
    #     outputfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed/kmeans_GSE169445_Nash_DC_combined"+str(index)+".txt"
    #     writesampletofilename(testallsample,outputfile)
    #     index+=1



    # samplenames=['AD1_AD2','AD3_AD4','AD5_AD6','Ct1_Ct2','Ct3_Ct4','Ct5_Ct6']
    # inputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/"
    # typefilepath= "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/scRNA_metadata.tsv"
    # samplecol=1
    # typecol = 8
    
    # for samplename in samplenames:
    #     allsamplefilepath = inputdir+samplename+".txt"
    #     outputdir=inputdir+samplename+"/"
    #     splitmatrixwithtype(allsamplefilepath,typefilepath,samplecol,typecol,outputdir,needtype=False)





    

    





