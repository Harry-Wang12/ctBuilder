
def filterlink(degenelist,graphdict):
    # print("Oringal graph has: ", str(len(graphdict)) ," item")
    filteredgraphdict={}
    for TF,targetlist in graphdict.items():
        if TF in degenelist:
            filteredgraphdict[TF]=[]
            for gene in targetlist:
                if gene in degenelist:
                    filteredgraphdict[TF].append(gene)
    # print("Filtered graph has: ", str(len(graphdict))  ," item")
    return filteredgraphdict



def loadDElist(degenelistfilepath):
    genelist= []
    with open(degenelistfilepath,"r") as degenelistfile:
        for line in degenelistfile.readlines():
            genelist.append(line.strip().upper())
    
    return genelist



def loadgraphfromfile(filepath,TFcol=1,target=2):
    
    graphdict = {}
    genelist = []
   

    with open(filepath,"r") as graphfile:
        for line in graphfile.readlines():
            linedata= line.strip().split("\t")
            tfgenename = linedata[TFcol-1]
            targetgenename = linedata[target-1]
            if not tfgenename in graphdict.keys():
                graphdict[tfgenename]=[]
                    
            if not targetgenename in graphdict[tfgenename]:
                graphdict[tfgenename].append(targetgenename)
            
            if not tfgenename in genelist:
                genelist.append(tfgenename)
            
            if not targetgenename in genelist:
                genelist.append(targetgenename)

    return  graphdict,genelist





