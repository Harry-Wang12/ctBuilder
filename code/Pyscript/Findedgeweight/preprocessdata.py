
from scipy.stats import ranksums


def findthedifferentgene(filepath1,filepath2,outputfilepath,threshold):
    content1={}
    content2={}
    genelist = []

    titleline = "Gene"
    # load file1

    with open(filepath1,'r') as file1:
        istitle=True

        for line in file1.readlines():
            linedata = line.strip().split("\t")
            if istitle:
                del(linedata[0])
                for titlename in linedata:
                    titleline+="\t"+titlename
                istitle=False
            else:
                
                genename = linedata[0].upper()
                del(linedata[0])
                content1[genename]=[]
                for genevalue in linedata:
                    content1[genename].append(float(genevalue))
                if not genename in genelist:
                    genelist.append(genename)
    
    with open(filepath2,'r') as file2:
        istitle=True

        for line in file2.readlines():
            linedata = line.strip().split("\t")
            if istitle:
                del(linedata[0])
                for titlename in linedata:
                    titleline+="\t"+titlename
                istitle=False
            else:
                
                genename = linedata[0].upper()
                del(linedata[0])
                content2[genename]=[]
                for genevalue in linedata:
                    content2[genename].append(float(genevalue))
                if not genename in genelist:
                    genelist.append(genename)
    with open(outputfilepath,'w') as outputfile:
        outputfile.write(titleline+'\n')
        for genename in genelist:
            if genename in content1.keys() and genename in content2.keys():
                staticvalue, pvalue  = ranksums(content1[genename], content2[genename])
                if pvalue<= threshold:
                    print(genename,pvalue)
                    pline=genename
                    for genevalue in content1[genename]:
                        pline+="\t"+str(genevalue)
                    for genevalue in content2[genename]:
                        pline+="\t"+str(genevalue)
                    
                    outputfile.write(pline+'\n')


    # with open(outputfile1path,'w') as outputfile:
    #     outputfile.write(titleline+'\n')
    #     for genename in genelist:
    #         if genename in content1.keys() and genename in content2.keys():
    #             staticvalue, pvalue  = ranksums(content1[genename], content2[genename])
    #             if pvalue<= threshold:
    #                 print(genename,pvalue)
    #                 pline=genename
    #                 for genevalue in content1[genename]:
    #                     pline+="\t"+str(genevalue)
    #                 # for genevalue in content2[genename]:
    #                 #     pline+="\t"+str(genevalue)
    #                 outputfile.write(pline+'\n')

    # with open(outputfile2path,'w') as outputfile:
    #     outputfile.write(titleline+'\n')
    #     for genename in genelist:
    #         if genename in content1.keys() and genename in content2.keys():
    #             staticvalue, pvalue  = ranksums(content1[genename], content2[genename])
    #             if pvalue<= threshold:
    #                 print(genename,pvalue)
    #                 pline=genename
    #                 # for genevalue in content1[genename]:
    #                 #     pline+="\t"+str(genevalue)
    #                 for genevalue in content2[genename]:
    #                     pline+="\t"+str(genevalue)
    #                 outputfile.write(pline+'\n')
    


if __name__=="__main__":

    filepath1="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/AD5_AD6/OPC.txt"
    filepath2="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/Ct1_Ct2/OPC.txt"
    threshold=0.1
    outputfilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/combinedranksum/0.1/GSE138852_AD_control_OPC_small_"+str(threshold)+".txt"
    findthedifferentgene(filepath1,filepath2,outputfilepath,threshold)



