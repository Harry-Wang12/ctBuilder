import numpy as np
import os
import statistics
import math
# this data is rescaled try to use the original first and then try to make the min become 0

# get split the cell and calculate log 2 without adding anything
def rescale(datalist,low,high):
    newdataset = []
    for value in datalist:
        newdataset.append((value - min(datalist)/(max(datalist)-min(datalist))*(high-low)))
    return newdataset

def unlog2(datalist):
    newdataset = []
    for value in datalist:
        newdataset.append(2**value)
    return newdataset


# rescale from 0 to 100
def generatelog2rescale(datatxtfile,outputdir,celltypefile,selectedcelltype,low,high):

    selectedcellname =[]

    with open(celltypefile,"r") as celltypes:
        title = True
        for line in celltypes.readlines():
            if title:
                title=False
            else:
                linedata=line.strip().split("\t")
                cellname = linedata[0]
                celltype = linedata[1]
                if celltype==selectedcelltype:
                    selectedcellname.append(cellname)
    
    samplelist={}
    allgene = []

    with open(datatxtfile,"r") as matrixfile:

            isheader = True
            filesamleplist={}
            headerorder=[]
            print("reading rescale")
            for line in matrixfile.readlines():
                if isheader:
                    linedata = line.strip().split("\t")
                    del(linedata[0])
                    for cellname in linedata:
                        filesamleplist[cellname]={}
                        headerorder.append(cellname)
                    isheader=False
                else:                 
                   
                    linedata = line.strip().split("\t")
                    genename = linedata[0].upper()
                    print(genename)
                    floatgenelist = []
                    del(linedata[0])
                    for i in range(len(linedata)):
                        floatgenelist.append(float(linedata[i]))

                    rescaleddlist = rescale(floatgenelist,low,high)
                    mean = sum(rescaleddlist)/len(rescaleddlist)
                    if not genename in allgene:
                        allgene.append(genename)                    
                    
                    for i in range(len(rescaleddlist)):
                        cellname = headerorder[i]
                        filesamleplist[cellname][genename]=np.log2((rescaleddlist[i])/mean)
                       
            for cellname,genelist in filesamleplist.items():
                print(cellname)
                if cellname in selectedcellname:
                    samplelist[cellname] = genelist
            del(filesamleplist)

    outputfile = outputdir+"combined.txt"
    with open(outputfile,"w") as output:
        print("writing rescale")
        title = "Sample"
        for cellname in samplelist.keys():
            title+="\t"+cellname
        output.write(title+"\n")
        for genename in allgene:
            pline = genename
            for cellname,genelist in samplelist.items():
                pline+='\t'+str(genelist[genename])
            output.write(pline+"\n")
    
    for cellname,genelist in samplelist.items():
        print(cellname)
        with open(outputdir+cellname.replace("-","_")+".txt","w") as indivalfile:
            title = "GENE\tRNA\n"
            indivalfile.write(title)
            for genename in allgene:
                pline = genename
                pline+='\t'+str(genelist[genename])
                indivalfile.write(pline+"\n")



        print(cellname)
        with open(outputdir+cellname.replace("-","_")+".txt","w") as indivalfile:
            title = "GENE\tRNA\n"
            indivalfile.write(title)
            for genename in allgene:
                pline = genename
                pline+='\t'+str(genelist[genename])
                indivalfile.write(pline+"\n")



def generatelog2testvsNormal(testdatatxtfile,normaldatatxtfile,outputfilepath):

    # read normal dataset 
    # allgene=[]

    normalaveragegenelist={}
    with open(normaldatatxtfile,"r") as matrixfile:
        isheader = True        
        print("reading normal")
        for line in matrixfile.readlines():
            if isheader:
               isheader=False
            else:
                linedata = line.strip().split("\t")
                genename = linedata[0].upper()
                
                # if not genename in allgene:
                #     allgene.append(genename)
                floatgenelist = []
                del(linedata[0])
                for i in range(len(linedata)):
                    floatgenelist.append(float(linedata[i]))
                # rescaleddlist = unlog2(floatgenelist)
                                
                # for i in range(len(rescaleddlist)):
                #     cellname = headerorder[i]
                #     filesamleplist[cellname][genename]=rescaleddlist[i]
                normalaveragegenelist[genename]=statistics.mean(floatgenelist)


    
    with open(outputfilepath,'w') as outputfile:
        with open(testdatatxtfile,"r") as matrixfile:            
            isheader = True
            print("reading test")
            for line in matrixfile.readlines():
                if isheader:
                    outputfile.write(line)
                    isheader=False
                else:
                    linedata = line.strip().split("\t")
                    genename = linedata[0].upper()
                    
                    
                    floatgenelist = []
                    del(linedata[0])

                    if genename in normalaveragegenelist.keys() and not "." in genename and not "-AS" in genename:
                        pline = genename
                        for genevalue in  linedata:
                            genevalue = float(genevalue)                            
                            logvalue = round(math.log2((genevalue+0.001)/(normalaveragegenelist[genename]+0.001)),2)
                            pline+='\t'+str(logvalue)
                        print(genename)
                        outputfile.write(pline+'\n')
                           

def splitbigmatrixtosplit(inputfilepath,outputdir):

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    samplelist ={}


    with open(inputfilepath,'r') as inputfile:
        isheader = True
        headerordr = []
        for line in inputfile.readlines():
            if isheader:
                linedata=line.strip().split('\t')
                del(linedata[0])
                headerordr=linedata
                for i in range(len(headerordr)):
                    outputfilepath = outputdir+headerordr[i]+'.txt'
                    samplelist[headerordr[i]]={}
                    with open(outputfilepath,'w') as outputfile:
                        outputfile.write("GENE\tRNA\n")

                isheader = False
            else:
                linedata=line.strip().split('\t')
                genename = linedata[0]
                print(genename)
                del(linedata[0])
                for i in range(len(headerordr)):
                    # 
                    samplelist[headerordr[i]][genename] = str(linedata[i])
                    # with open(outputfilepath,'a') as outputfile:
                    #     pline = genename+'\t'+str(linedata[i])+'\n'
                    #     outputfile.write(pline)

    print("writing")
    index = 0 
    for samplename, genenameandvalue in samplelist.items():
        index+=1
        outputfilepath = outputdir+samplename+'.txt'
        with open(outputfilepath,'a') as outputfile:
            print(outputfilepath,str(index/len(samplelist)))
            for genename,genevalue in genenameandvalue.items():
                pline = genename+'\t'+str(genevalue)+'\n'
                outputfile.write(pline)


                


                









# if __name__=="__main__":



#     # testdatatxtfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE5203285_NASH_DC_combined.txt"
#     # normaldatatxtfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSM5203286_control_DC_combined.txt"
#     outputfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE44770/LOAD/GSE44770_mapped.txt"
#     outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE44770/LOAD/indiv/"
#     # generatelog2testvsNormal(testdatatxtfile,normaldatatxtfile,outputfilepath)
#     splitbigmatrixtosplit(outputfilepath,outputdir)