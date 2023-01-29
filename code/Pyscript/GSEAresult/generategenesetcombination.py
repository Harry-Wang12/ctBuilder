

import kgmlparsar
import os
from itertools import permutations



def combinegenesets(geneset1,geneset2):
    returngenelist =[]
    for gene in geneset1:
        if not gene in returngenelist:
            returngenelist.append(gene)
    for gene in geneset2:
        if not gene in returngenelist:
            returngenelist.append(gene)
    
    return returngenelist


def generategmt(allgenedict,outputfilepath):
    
    with open(outputfilepath,'w') as outputfile:
        for description,genelist in allgenedict.items():
            pline = description+'\t'+description
            for gene in genelist:
                pline+='\t'+gene
            pline+='\n'
            outputfile.write(pline)

def readgenelistfromfile(filepath):
    genelist=[]
    with open(filepath, 'r') as file:
        for line in file.readlines():
            genename = line.strip()
            if not genename in genelist:
                genelist.append(genename)
    
    return genelist

def removesamegene(genelist1,genelist2):
    genelist = []

    for genevalue in genelist1:
        if not genevalue in genelist2:
            genelist.append(genevalue)
    return genelist







if __name__=="__main__":
    gene_id_symbol_map=kgmlparsar.getgeneidsymbolmap()

    dirpath= 'C:/Users/whl19/Documents/Code/GenebetweenPathways/kgmlfiles/'

    betweengenesetdir ="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_7-4/betweenroute_NASH/"

    graphlists={}
    genelists={}
    kgmlfilelist = []

    g = os.walk(dirpath)  

    for path,dir_list,file_list in g:  
        for file_name in file_list:  
            # kgmllist.append(os.path.join(path, file_name))
            pathwaykgmlfile=os.path.join(path, file_name)
            kgmlcontent = kgmlparsar.readkgml(pathwaykgmlfile,gene_id_symbol_map)
            graphdict,genelist = kgmlparsar.generateGraph(kgmlcontent,gene_id_symbol_map)
            kgmlfilelist.append(file_name)
            graphlists[file_name] = graphdict
            genelists[file_name] = genelist

    
    filecombinations = list(permutations(kgmlfilelist, 2))
    
    twopathwaysdict = {}

    for filecombination in filecombinations:
        description = filecombination[0]+'_'+filecombination[1]
        betweenfile = betweengenesetdir+filecombination[0]+'_'+filecombination[1]+'_combinelist.txt'
        betweengenelist = readgenelistfromfile(betweenfile)

        combinegenelist = combinegenesets(genelists[filecombination[0]],genelists[filecombination[1]])
        removegenelist = removesamegene(betweengenelist,combinegenelist)
        # combinegenelist = combinegenesets(combinegenelist,betweengenelist)
        if len(betweengenelist)>0:
            twopathwaysdict[description] = removegenelist
    

    
    # for filename,genelist in genelists.items():
    #     description = filename
    #     # betweenfile = betweengenesetdir+filecombination[0]+'_'+filecombination[1]+'_combinelist.txt'
    #     # betweengenelist = readgenelistfromfile(betweenfile)
    #     # combinegenelist = combinegenesets(genelists[filecombination[0]],genelists[filecombination[1]])
    #     # combinegenelist = combinegenesets(combinegenelist,betweengenelist)
    #     twopathwaysdict[description] = genelist


    outputfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/pathways/onlypathway_between7_8.gmt"

    generategmt(twopathwaysdict,outputfilepath)
