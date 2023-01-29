

# import Revisedpathway
# import pathwayprocess
# import dataprocess
# import CoexpressionDB
# import filteroutgene
import Dataprocess.GetRouteandCoexpression



if __name__ == "__main__":


    # # pathwayfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayexample.json"
    # # GRNfilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/GRNdb/SKCM_TCGA-regulons.txt"
    # # outputfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayexample_testrevise.txt"
    # # # Revisedpathway.revisepathway(pathwayfilepath,GRNfilepath,outputfile)
    # # Revisedpathway.reviseGRNhighpathway(pathwayfilepath,GRNfilepath,outputfile)
    # threshold=0.9
    # # pathwayfile1path = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/pathways/KEGG - TNF signaling pathway.json"
    # # pathwayfile2path = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/pathways/KEGG - MAPK signaling.json"
    # # GRNfilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/GRNdb/SKCM_TCGA-regulons.txt"
    # # outputfile1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Cohortedpathway/revised/TNF_GRN.txt"
    # # outputfile2 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Cohortedpathway/revised/MAPK_GRN.txt"

    # # outputfile3 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Cohortedpathway/revised/TNF_GRN_filtered.txt"
    # # outputfile4 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Cohortedpathway/revised/MAPK_GRN_filtered.txt"
    # # outputfile5 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Cohortedpathway/revised/MAPK_TNF_GRN_filtered_"+str(threshold)+".txt"
    # # outputfile6 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Cohortedpathway/revised/MAPK_TNF_KEGG_copexress.txt"
    # # outputfile7 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Cohortedpathway/revised/MAPK_TNF_KEGG.txt"
    # outputfile8 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Cohortedpathway/revised/MAPK_TNF_KEGG_IRF1_RELB.txt"
    # outputfile9 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Cohortedpathway/revised/MAPK_TNF_KEGG_IRF1_RELB_filter_"+str(threshold)+"_train.txt"



    # # # log2datafile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"

    # # # pathwayfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayexample.json"
    # # # biogridfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/BiogridDB/BIOGRID-ORGANISM-Homo_sapiens-4.2.192.tab3.txt"
    # # # Stringlinkfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Homo sampines/9606.protein.links.v11.0.txt"
    # # # Stringinfofilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Homo sampines/9606.protein.info.v11.0.txt"
   
    # cohortsupdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/IRF8up_down/uptrain/"

    # cohortsdowndir="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/IRF8up_down/downtrain/"

    # cohortsupdirtest="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/IRF8up_down/uptest/"

    # cohortsdowndirtest="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/IRF8up_down/downtest/"



    # # combinedpathway = Revisedpathway.combinepathways(pathwayfile1path,pathwayfile2path,pathwayname="Combined TNF and MAPK")

    # # log2data = dataprocess.loadfile(log2datafile,title=True)

    # # # Revisedpathway.reviseGRNhighpathway(pathwayfile1path,GRNfilepath,outputfile1)
    # # # Revisedpathway.reviseGRNhighpathway(pathwayfile2path,GRNfilepath,outputfile2)

    # # Revisedpathway.reviseGRNhighAndCorpathway(pathwayfile1path,GRNfilepath,log2data,outputfile3,threshold)
    # # Revisedpathway.reviseGRNhighAndCorpathway(pathwayfile2path,GRNfilepath,log2data,outputfile4,threshold)



    # # combinedpathway =  Revisedpathway.combinedtwoGRNpathway(outputfile3,outputfile4,pathwayname="Combined TNF MAPK and GRN")
    # # pathwayprocess.dumpjsontofile(outputfile7,combinedpathway)


    # # load biogrid
    

    
    # # Coexpresslist1 = CoexpressionDB.loadfileBiogrid(biogridfilepath)
    # # Coexpresslist2 = CoexpressionDB.loadfileString(Stringlinkfilepath,Stringinfofilepath)
    # # Coexpresslist = CoexpressionDB.combinetwodb(Coexpresslist1,Coexpresslist2)
    # # newpathwaydict = Revisedpathway.addtwoGenesCommonanduniquebundle("IRF1","RELB",outputfile7,Coexpresslist)

    # # # newpathwaydict = Revisedpathway.addbundletoTFCoexpress(outputfile7,Coexpresslist)
    # # # pathwayprocess.dumpjsontofile(outputfile6,newpathwaydict)
    # # # newpathwaydict = Revisedpathway.findcommongenebetweenBundles()

    # # pathwayprocess.dumpjsontofile(outputfile8,newpathwaydict)

    # # newpathwaydict = Revisedpathway.fliteroutBundlegene("IRF1_RELB",cohortsupdir,cohortsdowndir,outputfile8,threshold)

    # # pathwayprocess.dumpjsontofile(outputfile9,newpathwaydict)

    # # allgenelistfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/tmp/Combined TNF and MAPK_gene_names.txt"
    # # remainedgenelistfile="C:/Users/whl19/Documents/Code/GenebetweenPathways/tmp/Combined TNF and MAPK_gene_names_label.txt"
    # # traininggenelist = "C:/Users/whl19/Documents/Code/GenebetweenPathways/tmp/Combined TNF and MAPK_gene_names_small.txt"
    # # lr=0.0001
    # # epchos=10000
    # # nnoutputfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/achieve/result5.txt"
    # # filteroutgene.Method1(cohortsupdir,cohortsdowndir,cohortsupdirtest,cohortsdowndirtest,allgenelistfile,traininggenelist,remainedgenelistfile,lr,epchos,nnoutputfile,TFgene="IRF1")



    
    # # filteroutgene.randomdatasplit(10,"C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_raw.txt","C:/Users/whl19/Documents/Code/GenebetweenPathways/tmp/TGFA_CCL2_gene_names.txt","C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/randomselect/test7_ccl2/")
    # # # filteroutgene.usingGENIN3getresult("C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/randomselect/test1/",8)

    # # findname="CCL2 

    # # dataprocess.tranposmatrix("C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt","C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio_T.txt")
    # # dataprocess.tranposmatrix("C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_raw.txt","C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_raw_T.txt")













    
    inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    outputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-18-2021/TCGA_SKCM_RouteScore.txt"
    CoexpressionDir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-18-2021/CoexpressionDir/"

    Dataprocess.GetRouteandCoexpression.generateRouteAndCoexpressionfileDown(inputtable, outputtable,CoexpressionDir)
    