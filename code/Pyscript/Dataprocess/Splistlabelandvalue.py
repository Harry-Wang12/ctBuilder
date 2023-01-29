import random

def loadfile(allsamplefilepath):
    
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