# Split TCGA data with its cancer stage
# And some delete gene with std lower than 1.



import statistics



if __name__ == "__main__":
    Expressiondatapath="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_raw.txt"
    Clincialdatapath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_samples.txt"
    outputfolder = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/filteredStage/"
    MINaverage=5
    MINstdev = 5


    SampleMap={}
    expressiondata = {}
    expressiontitle = []
    genelist = []
    

    with open(Clincialdatapath,"r") as Clincialdata:
        title = True
        for line in Clincialdata.readlines():
            if title:
                title = False
            else:
                linedata = line.strip().split("\t")
                Samplename = linedata[0]
                Stage=linedata[3].replace("/","_")
                if Stage not in SampleMap.keys():
                    SampleMap[Stage] = []
                    SampleMap[Stage].append(Samplename)
                else:
                    SampleMap[Stage].append(Samplename)

    
    with open(Expressiondatapath,"r") as Expressiondata:
        title = True
        for line in Expressiondata.readlines():
            if title:
                linedata = line.strip().split("\t")
                del(linedata[0])
                expressiontitle = linedata
                for titlename in expressiontitle:
                    expressiondata[titlename]=[]
                title = False
            else:
                linedata = line.strip().split("\t")
                genename = linedata[0]
                
                del(linedata[0])
                linestdev = statistics.stdev([float(i) for i in linedata])
                lineavg = statistics.mean([float(i) for i in linedata])
                if lineavg > MINaverage and linestdev > MINstdev:
                        # print(genename)
                    genelist.append(genename)
                    for i in range(len(linedata)):
                        expressiondata[expressiontitle[i]].append(linedata[i])
    

    for key,titlenames in SampleMap.items():
        print(key)
        # print(titlenames)
        with open(outputfolder+key+"_"+str(MINaverage)+"_"+str(MINstdev)+".txt","w") as outputfile:
            ptitle = "samplename"
            for genename in genelist:
                ptitle+="\t"+genename
            ptitle+="\n"
            outputfile.write(ptitle)
            for title in titlenames:
                pline = title
                for value in expressiondata[title]:
                    pline+="\t"+value
                pline+="\n"
                outputfile.write(pline)




                
                






