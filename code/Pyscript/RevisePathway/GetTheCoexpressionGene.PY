
import CoexpressionDB




def getCoexpressionGene(genename):
    biogridfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/BiogridDB/BIOGRID-ORGANISM-Homo_sapiens-4.2.192.tab3.txt"
    Stringlinkfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Homo sampines/9606.protein.links.v11.0.txt"
    Stringinfofilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Homo sampines/9606.protein.info.v11.0.txt"
    Coexpresslist1 = CoexpressionDB.loadfileBiogrid(biogridfilepath)
    Coexpresslist2 = CoexpressionDB.loadfileString(Stringlinkfilepath,Stringinfofilepath)
    Coexpresslist = CoexpressionDB.combinetwodb(Coexpresslist1,Coexpresslist2)
    Coexpressionlist= []

    for genepair in Coexpresslist:
        if genepair["gene1"] == genename and not genepair["gene2"] in Coexpressionlist:
            Coexpressionlist.append(genepair["gene2"])
        elif genepair["gene2"] == genename and not genepair["gene1"] in Coexpressionlist:
            Coexpressionlist.append(genepair["gene1"])
    
    return Coexpressionlist
            










