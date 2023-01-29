
import random
import os

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

def seperateall(numberofgroup,allsamples,israndom):

    if israndom:
        random.shuffle(allsamples)

    seperatedsamples=[]

    for i in range(numberofgroup):
        seperatedsamples.append([])

    index =-1
    for sample in allsamples:
        if index==numberofgroup-1:
            index =0
        else:
            index+=1
            
        seperatedsamples[index].append(sample)      

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






def generatetmpfile(allsamples,outputfile):
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

    samplenames=['AD1_AD2','AD3_AD4','AD5_AD6','Ct1_Ct2','Ct3_Ct4','Ct5_Ct6']
    inputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/"
    typefilepath= "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/scRNA_metadata.tsv"
    samplecol=1
    typecol = 8
    
    for samplename in samplenames:
        allsamplefilepath = inputdir+samplename+".txt"
        outputdir=inputdir+samplename+"/"
        splitmatrixwithtype(allsamplefilepath,typefilepath,samplecol,typecol,outputdir,needtype=False)





    

    





