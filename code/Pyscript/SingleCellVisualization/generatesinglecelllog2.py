import numpy as np
import os
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


def generatelog2unlog(datatxtfile,outputdir,celltypefile,selectedcelltype):
    outputfile = outputdir+"combined.txt"
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
            print("reading unlog")
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

                    rescaleddlist = unlog2(floatgenelist)
                    mean = sum(rescaleddlist)/len(rescaleddlist)
                    if not genename in allgene:
                        allgene.append(genename)
                    
                    
                    for i in range(len(rescaleddlist)):
                        cellname = headerorder[i]
                        filesamleplist[cellname][genename]=np.log2((rescaleddlist[i])/mean)
                       
            for cellname,genelist in filesamleplist.items():
                if cellname in selectedcellname:
                    samplelist[cellname] = genelist
            del(filesamleplist)

    with open(outputfile,"w") as output:
        print("writing unlog")
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



def generatelog2unlogNashvsNormal(nashdatatxtfile,normaldatatxtfile,outputfiledir,nashcelltypefile,normalcelltypefile,selectedcelltype):

    normalselectedcellname =[]

    with open(normalcelltypefile,"r") as celltypes:
        title = True
        for line in celltypes.readlines():
            if title:
                title=False
            else:
                linedata=line.strip().split("\t")
                cellname = linedata[0]
                celltype = linedata[1]
                if celltype==selectedcelltype:
                    normalselectedcellname.append(cellname)
    allgene=[]
    normalaveragegenelist={}
    with open(normaldatatxtfile,"r") as matrixfile:
        isheader = True
        filesamleplist={}
        headerorder=[] 
        print("reading normal")
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
                rescaleddlist = unlog2(floatgenelist)
                
                if not genename in allgene:
                    allgene.append(genename)
                for i in range(len(rescaleddlist)):
                    cellname = headerorder[i]
                    filesamleplist[cellname][genename]=rescaleddlist[i]

        for gene in allgene:
            normalaveragegenelist[gene]=[]          
            for cellname,genelist in filesamleplist.items():
                if cellname in normalselectedcellname:
                    normalaveragegenelist[gene].append(genelist[genename])
        for gene in allgene:
            normalaveragegenelist[gene]=sum(normalaveragegenelist[gene])/len(normalaveragegenelist[gene])

        del(filesamleplist)

    nashlog2samplelist={}
    nashallgene=[]
    nashselectedcellname =[]
    with open(nashcelltypefile,"r") as celltypes:
        title = True
        for line in celltypes.readlines():
            if title:
                title=False
            else:
                linedata=line.strip().split("\t")
                cellname = linedata[0]
                celltype = linedata[1]
                if celltype==selectedcelltype:
                    nashselectedcellname.append(cellname)

    with open(nashdatatxtfile,"r") as matrixfile:            
        isheader = True
        filesamleplist={}
        headerorder=[] 
        print("reading unlog")
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

                rescaleddlist = unlog2(floatgenelist)
                if genename in normalaveragegenelist.keys():
                    mean = normalaveragegenelist[genename]  

                    if not genename in nashallgene:
                        nashallgene.append(genename)       
                    
                    for i in range(len(rescaleddlist)):
                        cellname = headerorder[i]
                        filesamleplist[cellname][genename]=np.log2((rescaleddlist[i])/mean)
                        
        for cellname,genelist in filesamleplist.items():
            if cellname in nashselectedcellname:
                nashlog2samplelist[cellname] = genelist
        del(filesamleplist)

    outputfilesummary=outputfiledir+"combined.txt"
    with open(outputfilesummary,"w") as output:
        print("writing unlog")
        title = "Sample"
        for cellname in nashlog2samplelist.keys():
            title+="\t"+cellname
        output.write(title+"\n")

        for genename in nashallgene:
            if genename in normalaveragegenelist.keys():
                pline = genename
                for cellname,genelist in nashlog2samplelist.items():
                    pline+='\t'+str(genelist[genename])
                output.write(pline+"\n")

    for cellname,genelist in nashlog2samplelist.items():
        print(cellname)
        with open(outputfiledir+cellname.replace("-","_")+".txt","w") as indivalfile:
            title = "GENE\tRNA\n"
            indivalfile.write(title)
            for genename in nashallgene:
                pline = genename
                pline+='\t'+str(genelist[genename])
                indivalfile.write(pline+"\n")



def generatelog2unlogNashvsNormalforGSEA(nashdatatxtfile,normaldatatxtfile,outputfilepath,nashcelltypefile,clsfilepath,normalcelltypefile,selectedcelltype):
    normalselectedcellname =[]
    with open(normalcelltypefile,"r") as celltypes:
        title = True
        for line in celltypes.readlines():
            if title:
                title=False
            else:
                linedata=line.strip().split("\t")
                cellname = linedata[0]
                celltype = linedata[1]
                if celltype==selectedcelltype:
                    normalselectedcellname.append(cellname)
    normalallgene=[]
    normalsamplelist={}
    
    with open(normaldatatxtfile,"r") as matrixfile:
        isheader = True
        filesamleplist={}
        headerorder=[] 
        print("reading normal")
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
                rescaleddlist = unlog2(floatgenelist)
                
                if not genename in normalallgene:
                    normalallgene.append(genename)
                for i in range(len(rescaleddlist)):
                    cellname = headerorder[i]
                    filesamleplist[cellname][genename]=rescaleddlist[i]

        for cellname,genelist in filesamleplist.items():
            if cellname in normalselectedcellname:
                normalsamplelist[cellname] = genelist
        

        del(filesamleplist)

    nashsamplelist={}
    nashallgene=[]
    nashselectedcellname =[]
    with open(nashcelltypefile,"r") as celltypes:
        title = True
        for line in celltypes.readlines():
            if title:
                title=False
            else:
                linedata=line.strip().split("\t")
                cellname = linedata[0]
                celltype = linedata[1]
                if celltype==selectedcelltype:
                    nashselectedcellname.append(cellname)

    with open(nashdatatxtfile,"r") as matrixfile:            
        isheader = True
        filesamleplist={}
        headerorder=[] 
        print("reading unlog")
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

                rescaleddlist = unlog2(floatgenelist)              

                if not genename in nashallgene:
                    nashallgene.append(genename)       
                
                for i in range(len(rescaleddlist)):
                    cellname = headerorder[i]
                    filesamleplist[cellname][genename]=rescaleddlist[i]
                        
        for cellname,genelist in filesamleplist.items():
            if cellname in nashselectedcellname:
                nashsamplelist[cellname] = genelist
        del(filesamleplist)


    combinegene = [ i for i in nashallgene if i in normalallgene]

    combinedsamplelist = {}

    with open(clsfilepath,"w") as clsfile:
        pline = str(len(normalsamplelist)+len(nashsamplelist))+" 2 1\n"
        pline+="# normal nash\n"
        for i in range(len(normalsamplelist)):
            pline+="0 "
        for i in range(len(nashsamplelist)):
            pline+="1 "
        
        clsfile.write(pline[:-1])


   
    for cellname,genelist in normalsamplelist.items():
        combinedsamplelist[cellname]={}
        for genename in combinegene:
            combinedsamplelist[cellname][genename] = genelist[genename]
    
    for cellname,genelist in nashsamplelist.items():
        combinedsamplelist[cellname]={}
        for genename in combinegene:
            combinedsamplelist[cellname][genename] = genelist[genename]
                
   
    with open(outputfilepath,"w") as output:
        print("writing unlog")
        title = "Sample"
        for cellname in combinedsamplelist.keys():
            title+="\t"+cellname
        output.write(title+"\n")

        for genename in combinegene:
            
            pline = genename
            for cellname,genelist in combinedsamplelist.items():
                pline+='\t'+str(genelist[genename])
            output.write(pline+"\n")






def generatelog2raw(datatxtfile,outputdir):

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)


    outputfile = outputdir+"combined.txt"
    # selectedcellname =[]

    # with open(celltypefile,"r") as celltypes:
    #     title = True
    #     for line in celltypes.readlines():
    #         if title:

    #             title=False
    #         else:
    #             linedata=line.strip().split("\t")
    #             cellname = linedata[0]
    #             if len(linedata)>11:
    #                 celltype = linedata[11]
    #             else:
    #                 celltype = "unknown"
    #             if celltype==selectedcelltype:
    #                 selectedcellname.append(cellname)
    
    samplelist={}
    allgene = []

    with open(datatxtfile,"r") as matrixfile:            
            isheader = True
            filesamleplist={}
            headerorder=[] 
            print("reading")
            for line in matrixfile.readlines():
                if isheader:
                    linedata = line.strip().split("\t")
                    # del(linedata[0])
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

                    rescaleddlist = floatgenelist
                    mean = sum(rescaleddlist)/len(rescaleddlist)
                    if not genename in allgene:
                        allgene.append(genename)                  
                    
                    for i in range(len(rescaleddlist)):
                        cellname = headerorder[i]
                        filesamleplist[cellname][genename]=np.log2((rescaleddlist[i]+0.01)/(mean+0.01))
                       
            for cellname,genelist in filesamleplist.items():
                # if cellname in selectedcellname:
                samplelist[cellname] = genelist
            del(filesamleplist)

    with open(outputfile,"w") as output:
        print("writing")
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






def generatelog2unlogNashvsNormal(datmatrixfile,outputfiledir,testcelltype,normalcelltype,celltypefile):

    normalselectedcellname =[]
    testselectedcellname = []

    with open(celltypefile,"r") as celltypes:
        title = True
        for line in celltypes.readlines():
            if title:
                title=False
            else:
                linedata=line.strip().split("\t")
                cellname = linedata[0]
                celltype = linedata[1]
                if celltype==normalcelltype:
                    normalselectedcellname.append(cellname)
                if celltype==testcelltype:
                    testselectedcellname.append(cellname)
    allgene=[]
    
    nashlog2samplelist={}
    with open(datmatrixfile,"r") as matrixfile:
        isheader = True
        filesamleplist={}
        headerorder=[] 
        print("reading normal")
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
                rescaleddlist = unlog2(floatgenelist)
                
                if not genename in allgene:
                    allgene.append(genename)
                for i in range(len(rescaleddlist)):
                    cellname = headerorder[i]
                    filesamleplist[cellname][genename]=rescaleddlist[i]

    normalaveragegenelist={}
    # calculate normal average
    for genename in allgene:
        genealllist =[]
        for cellname,genevaluelist in filesamleplist.items:
            if cellname in normalselectedcellname:
                genealllist.append(genevaluelist[genename])
        normalaveragegenelist[genename] = sum(genealllist)/len(genealllist)
           


    





        # for gene in allgene:
        #     normalaveragegenelist[gene]=[]          
        #     for cellname,genelist in filesamleplist.items():
        #         if cellname in normalselectedcellname:
        #             normalaveragegenelist[gene].append(genelist[genename])
        # for gene in allgene:
        #     normalaveragegenelist[gene]=sum(normalaveragegenelist[gene])/len(normalaveragegenelist[gene])

        # del(filesamleplist)

    # nashlog2samplelist={}
    # nashallgene=[]
    # nashselectedcellname =[]
    # with open(nashcelltypefile,"r") as celltypes:
    #     title = True
    #     for line in celltypes.readlines():
    #         if title:
    #             title=False
    #         else:
    #             linedata=line.strip().split("\t")
    #             cellname = linedata[0]
    #             celltype = linedata[1]
    #             if celltype==selectedcelltype:
    #                 nashselectedcellname.append(cellname)

    # with open(nashdatatxtfile,"r") as matrixfile:            
    #     isheader = True
    #     filesamleplist={}
    #     headerorder=[] 
    #     print("reading unlog")
    #     for line in matrixfile.readlines():
    #         if isheader:
    #             linedata = line.strip().split("\t")
    #             del(linedata[0])
    #             for cellname in linedata:
    #                 filesamleplist[cellname]={}
    #                 headerorder.append(cellname)
    #             isheader=False
    #         else:

    #             linedata = line.strip().split("\t")
    #             genename = linedata[0].upper()
    #             print(genename)
    #             floatgenelist = []
    #             del(linedata[0])
    #             for i in range(len(linedata)):
    #                 floatgenelist.append(float(linedata[i]))

    #             rescaleddlist = unlog2(floatgenelist)
    #             if genename in normalaveragegenelist.keys():
    #                 mean = normalaveragegenelist[genename]  

    #                 if not genename in nashallgene:
    #                     nashallgene.append(genename)       
                    
    #                 for i in range(len(rescaleddlist)):
    #                     cellname = headerorder[i]
    #                     filesamleplist[cellname][genename]=np.log2((rescaleddlist[i])/mean)
                        
    #     for cellname,genelist in filesamleplist.items():
    #         if cellname in nashselectedcellname:
    #             nashlog2samplelist[cellname] = genelist
    #     del(filesamleplist)

    # outputfilesummary=outputfiledir+"combined.txt"
    # with open(outputfilesummary,"w") as output:
    #     print("writing unlog")
    #     title = "Sample"
    #     for cellname in nashlog2samplelist.keys():
    #         title+="\t"+cellname
    #     output.write(title+"\n")

    #     for genename in nashallgene:
    #         if genename in normalaveragegenelist.keys():
    #             pline = genename
    #             for cellname,genelist in nashlog2samplelist.items():
    #                 pline+='\t'+str(genelist[genename])
    #             output.write(pline+"\n")

    # for cellname,genelist in nashlog2samplelist.items():
    #     print(cellname)
    #     with open(outputfiledir+cellname.replace("-","_")+".txt","w") as indivalfile:
    #         title = "GENE\tRNA\n"
    #         indivalfile.write(title)
    #         for genename in nashallgene:
    #             pline = genename
    #             pline+='\t'+str(genelist[genename])
    #             indivalfile.write(pline+"\n")







if __name__=="__main__":




    # nash
    # nashdatatxtfile =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/GSE155182_NASH_7m_liver.txt"
    # # outputfilerescale =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/GSE155182_Normal_7m_liver_rescale.txt"
    # NASHoutputdirunlog2 =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/NASH_POPULATION_unlog/"
    # nashcelltypefile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/GSE155182_NASH_7m_liver_label.txt"




    # normal
    # normaldatatxtfile =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/GSE155182_Normal_7m_liver.txt"
    # outputdirrescaleNASH =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/NASH_POPULATION_rescale/"
    # outputdirrescaleNormal =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/Normal_POPULATION_rescale/"
    # normaloutputdirunlog2 =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/Normal_POPULATION_unlog/"
    # normalcelltypefile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/GSE155182_Normal_7m_liver_label.txt"


    # outputfiledir =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/inflammatory_non-inflammatory_macrophage/"
    # selectedcelltype = "Macrophage"
    # low=0
    # high = 10
    # print("rescale")
    # generatelog2rescale(nashdatatxtfile,outputdirrescaleNASH,nashcelltypefile,selectedcelltype,low,high)
    # generatelog2rescale(normaldatatxtfile,outputdirrescaleNormal,normalcelltypefile,selectedcelltype,low,high)

    # generatelog2unlog(normaldatatxtfile,normaloutputdirunlog2,normalcelltypefile,selectedcelltype)
    # generatelog2unlog(nashdatatxtfile,NASHoutputdirunlog2,nashcelltypefile,selectedcelltype)



    # print("unlog")
    # generatelog2unlog(datatxtfile,outputfileunlog2,celltypefile,selectedcelltype)

    generatelog2unlogNashvsNormal(nashdatatxtfile,normaldatatxtfile,outputfiledir,nashcelltypefile,normalcelltypefile,selectedcelltype)

    # outputfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/unlog_NASH_NORMAL_GSEA.txt"
    # clsfilepath =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/unlog_NASH_NORMAL_GSEA.cls"
    # generatelog2unlogNashvsNormalforGSEA(nashdatatxtfile,normaldatatxtfile,outputfilepath,nashcelltypefile,clsfilepath,normalcelltypefile,selectedcelltype)










    # datatxtfile='C:/Users/whl19/Documents/Code/GenebetweenPathways/GSE115469/matrix_mapped.txt'
    # celltypefile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/GSE115469/cellclabel.txt'
    # selectedcelltype="inflammatory macrophage"
    # outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/GSE115469/"+selectedcelltype.replace(" ","_")+"/"
    # # generatelog2raw(datatxtfile,outputdir,celltypefile,selectedcelltype)


    # selectedcelltype="non-inflammatory macrophage"
    # outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/GSE115469/"+selectedcelltype.replace(" ","_")+"/"
    # generatelog2raw(datatxtfile,outputdir,celltypefile,selectedcelltype)


    for i in range(20):
        outdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/singlecellsimulation/"
        datatxtfile=outdir+"counts_"+str(i)+".txt"
        outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/singlecellsimulation/counts_"+str(i)+"_log2/"
        generatelog2raw(datatxtfile,outputdir)
