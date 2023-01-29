


import GetRouteandCoexpression
import linkpathway
import os





if __name__ == '__main__':

    
    importancethreshold = 0.05
    subgouppercentage = 0.8
    subgroupnumber = 10
    # 3-15 inflammatory_macrophage non_inflammatory_macrophage
    allresultdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/"
    inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/non-inflammatory_macrophage/combined.txt"
    # allresultdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-15-2021_GSE115469_non-inflamtory/"
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/non-inflammatory_macrophage/combined.txt"


    outputtable = allresultdir+"RouteScore.txt"
    CoexpressionDir = allresultdir+"CoexpressionDir/"
    subgroupdir = allresultdir+"CoexpressionSubgroup/"
    resultdir = allresultdir+"Revisedpathway/"  
    summaryfilepath =  subgroupdir+"summary.txt"
    pathwaygenedir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/pathways/Keggpathwaygene/"    
    pathwaydir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/pathways/"
    revisedScoring=allresultdir+"revisedScoring/"
   

    if not os.path.exists(allresultdir):
        os.makedirs(allresultdir)
    if not os.path.exists(CoexpressionDir):
        os.makedirs(CoexpressionDir)
    if not os.path.exists(subgroupdir):
        os.makedirs(subgroupdir)
    if not os.path.exists(resultdir):
        os.makedirs(resultdir)
    if not os.path.exists(revisedScoring):
        os.makedirs(revisedScoring)    


    GetRouteandCoexpression.generateRouteAndCoexpressionfileDown(inputtable, outputtable,CoexpressionDir,subgroupdir,subgroupnumber,subgouppercentage,importancethreshold)


    linkpathway.getsummaryofsubgroup(subgroupdir)
    
    # findlinkegene(summaryfilepath,pathwaygenedir,resultdir)

    linkpathway.findlinkpathwayforDownStreamdata(summaryfilepath,pathwaygenedir,pathwaydir,resultdir)





    








    # os.system("shutdown -s -t  1")