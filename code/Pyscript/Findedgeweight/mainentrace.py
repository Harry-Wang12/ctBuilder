
import geneexpressionprocess
import kgmlparsar
import featureselectiongiveweight
import statistics
import multiprocessing
import os
import postanalysis
import preprocessdata
import copy



if  __name__=="__main__":

    
    # init
    # allsamplefilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/matrix_mapped.txt"
    # typefilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/cellclabel.txt"
    # outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/"
    # tmpdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/tmp/"
    
    # cancerlog2data = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    

    # inflamtionlog2file = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/inflammatory_macrophage/combined.txt"

    # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allcontent_new.txt'
    
    gene_id_symbol_map = kgmlparsar.getgeneidsymbolmap()
    
    # celltype = "MCDD_Liver_c57bl_cDC_2wkfile1"

    # samplecol = 1
    # typecol = 11

    # importancethreshold = -1


    # print("generate expression file")

    # geneexpressionprocess.splitmatrixwithtype(allsamplefilepath,typefilepath,samplecol,typecol,outputdir,celltype)



    # print("loading expression file")

    # # load expression file

    # # allsamples = geneexpressionprocess.loadtotalfile(outputdir+celltype+".txt")

    # allsamples = geneexpressionprocess.loadtotalfile(MCDD_Liver_c57bl_cDC_2wkfile1)



    ########################################################################################
    # actually file
    #     
    ########################################################################################

    

                        



    # resultending ="MCDD_Liver_c57bl_cDC_2wkfile"
    # resultweightfile1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/result1_"+resultending+".txt"
    # resultweightfile2 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/result2_"+resultending+".txt"
    # MCDD_Liver_c57bl_cDC_2wksdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/MCDD_Liver_c57bl_cDC_2wks/"
    # MCDD_Liver_c57bl_cDC_2wkfile1 = MCDD_Liver_c57bl_cDC_2wksdir+"GSM5203271_AB2641.txt"
    # MCDD_Liver_c57bl_cDC_2wkfile2 = MCDD_Liver_c57bl_cDC_2wksdir+"GSM5203272_AB2642.txt"


    # # split into test1 and test2
    # print("loading expression file")

    # # load expression file

    # # allsamples = geneexpressionprocess.loadtotalfile(outputdir+celltype+".txt")

    # testsample1 = geneexpressionprocess.loadtotalfile(MCDD_Liver_c57bl_cDC_2wkfile1)
    # testsample2 = geneexpressionprocess.loadtotalfile(MCDD_Liver_c57bl_cDC_2wkfile2)
    # israndom = True
    # numberofgroup = 2
    # # Alltestsamples = geneexpressionprocess.seperateall(numberofgroup,allsamples,israndom)
    # # testsample1 = Alltestsamples[0]
    # # testsample2 = Alltestsamples[1]


    

    # print("loading graph")
    # # start load graph

    # # kgmldirpath= 'C:/Users/whl19/Documents/Code/GenebetweenPathways/kgmlfiles/'
    # # allkgmlcontent = kgmlparsar.loadkgmlfromdir(kgmldirpath)

    # # StringDBfile ="C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Homo sampines/9606.protein.links.v11.0.txt"
    # # Stringinfofile ="C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/StringDB/Homo sampines/9606.protein.info.v11.0.txt"
    # # BioGridfile="C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/BiogridDB/BIOGRID-ORGANISM-Homo_sapiens-4.2.192.tab3.txt"

    # # PPIcontent = kgmlparsar.loadPPIcontent(StringDBfile,Stringinfofile,BioGridfile,gene_id_symbol_map)



    # # allcontent = kgmlparsar.combinePPIdatabase(allkgmlcontent,PPIcontent)

    # # allkgmlcontent = kgmlparsar.combineContent(allkgmlcontent)

    # # allkgmlcontent = kgmlparsar.loadjsontocontent(jsonfile)
    # # sourcelist = kgmlparsar.findsourceofgene(allkgmlcontent)
    # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allsourcelink_new.txt'
    # sourcelist = kgmlparsar.loadjsontocontent(jsonfile)
    
    # # 

    # # start feature selection
    # iter = 1
    


    # subgroupnumber = 5

    
    # # # testsample1
    # # subgroupsample1 = geneexpressionprocess.seperateall(subgroupnumber,testsample1,israndom)
    # # for i in range(subgroupnumber):
    # #     tmpfile = tmpdir+"test1_"+str(i)+".txt"
    # #     geneexpressionprocess.generatetmpfile(subgroupsample1[i],tmpfile)
    
    # # # testsample2
    # # subgroupsample2 = geneexpressionprocess.seperateall(subgroupnumber,testsample2,israndom)
    # # for i in range(subgroupnumber):
    # #     tmpfile = tmpdir+"test2_"+str(i)+".txt"
    # #     geneexpressionprocess.generatetmpfile(subgroupsample2[i],tmpfile)

    # # manager = multiprocessing.Manager()
    # # weight1 = manager.dict()
    # # weight2 = manager.dict()
    # # threads =16
    # weight1={}
    # weight2={}

    # for i in range(iter):

    #     print("iteration "+str(i))

    #     # sperate
    #     subgroupsample1 = geneexpressionprocess.seperateall(subgroupnumber,testsample1,israndom)
    #     subgroupsample2 = geneexpressionprocess.seperateall(subgroupnumber,testsample2,israndom)

    #     geneindex= 0
    #     for targetgeneid,featuregeneids in sourcelist.items(): 
                    
    #         if targetgeneid in gene_id_symbol_map.keys():
    #             targetgene = gene_id_symbol_map[targetgeneid]
    #             geneindex+=1  
    #             # if geneindex%100 ==0:
    #             print(("gene id "+str(geneindex))+" "+targetgene +" : "+str(len(sourcelist)))
    #             if not targetgene  in weight1.keys():
    #                 weight1[targetgene]={}
    #             if not targetgene  in weight2.keys():
    #                 weight2[targetgene]={}
    #             featuregenelist = []
    #             for nodeid in featuregeneids:
    #                 if nodeid in gene_id_symbol_map.keys():
    #                     featuregenelist.append(gene_id_symbol_map[nodeid])            
    #             # for testgroup1
    #             # index = 0
    #             for samples in subgroupsample1:
    #                 # index+=1
    #                 # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
    #                 weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #                 for featurename,weight in weightlist.items():
    #                     if not featurename in weight1[targetgene].keys():
    #                         weight1[targetgene][featurename]=[weight]
    #                     else:
    #                         weight1[targetgene][featurename].append(weight)
    #             # index = 0
    #             for samples in subgroupsample2:
    #                 # index+=1
    #                 # print("subgroupnumber2 "+str(index)+":"+str(subgroupnumber))
    #                 weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #                 for featurename,weight in weightlist.items():
    #                     if not featurename in weight2[targetgene].keys():
    #                         weight2[targetgene][featurename]=[weight]
    #                     else:
    #                         weight2[targetgene][featurename].append(weight)

    # # after process
    # # need to figure out later, right now just calculate mean of each edge
    # with open(resultweightfile1,'w') as resultweightfile:
    #     for target,featureweights in weight1.items():
    #         for featurename,featureweight in featureweights.items():
    #             if statistics.mean(featureweight) >= importancethreshold:
    #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                 resultweightfile.write(pline)
    
    # with open(resultweightfile2,'w') as resultweightfile:
    #     for target,featureweights in weight2.items():
    #         for featurename,featureweight in featureweights.items():
    #             if statistics.mean(featureweight) >= importancethreshold:
    #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                 resultweightfile.write(pline)
            
########################################################################################
    # actually file
    #     
########################################################################################

########################################################################################
# Generate random file
#     
########################################################################################
    # randomfile='C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/singlecellsimulation/new_counts_0.txt'
    # print("loading expression file")

    # # load expression file

  
    # allsamples = geneexpressionprocess.loadtotalfile(randomfile)

    # israndom = True
    # subgroupnumber = 50
    # iter = 1
    # print("loading source list")
    # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allsourcelink_new_allgene_test.txt'
    # sourcelist = kgmlparsar.loadjsontocontent(jsonfile)
    

    # weight1={}
    
    # for i in range(iter):
    #     print("iteration "+str(i))
    #     # sperate
    #     subgroupsample1 = geneexpressionprocess.seperateall(subgroupnumber,allsamples,israndom)
    #     geneindex= 0
    #     for targetgeneid,featuregeneids in sourcelist.items():                         
    #         if targetgeneid in gene_id_symbol_map.keys():
    #             targetgene = gene_id_symbol_map[targetgeneid]
    #             geneindex+=1  
    #             # if geneindex%100 ==0:
    #             print(("gene id "+str(geneindex))+" "+targetgene +" : "+str(len(sourcelist)))
    #             if not targetgene  in weight1.keys():
    #                 weight1[targetgene]={}
    #             featuregenelist = []
    #             for nodeid in featuregeneids:
    #                 if nodeid in gene_id_symbol_map.keys():
    #                     featuregenelist.append(gene_id_symbol_map[nodeid])            
    #             # for testgroup1
    #             # index = 0
    #             for samples in subgroupsample1:
    #                 # index+=1
    #                 # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
    #                 weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #                 for featurename,weight in weightlist.items():
    #                     if not featurename in weight1[targetgene].keys():
    #                         weight1[targetgene][featurename]=[weight]
    #                     else:
    #                         weight1[targetgene][featurename].append(weight)
    
    # randomresultfile1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/new_result_pvalue_withppi_allgene_test.txt"
    # with open(randomresultfile1,'w') as resultweightfile:
    #     for target,featureweights in weight1.items():
    #         for featurename,featureweight in featureweights.items():                
    #                 pline = target+"\t"+featurename
    #                 for weight in featureweight:
    #                     pline+="\t"+str(weight)
    #                 resultweightfile.write(pline+'\n')            

    # os.system("shutdown -s -t  1")

    ########################################################################################
    # Generate random file
    #     
    ########################################################################################

    ########################################################################################
    # Post analysis
    #     
    ########################################################################################
    # resultending ="MCDD_Liver_c57bl_cDC_2wkfile"
    # resultweightfile1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/result1_"+resultending+".txt"
    # resultweightfile2 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/result2_"+resultending+".txt"
    # # randomfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/result_pvalue_withppi.txt"
    # # randomnormaldic = postanalysis.loadrandomfile(randomfilepath)

    # # postanalysis.reviseresultfile(resultweightfile1,randomnormaldic)
    # # postanalysis.reviseresultfile(resultweightfile2,randomnormaldic)
    
    # outputfilepath1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/testresult/result1_"+resultending+"_refine.txt"
    # outputfilepath2 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/testresult/result2_"+resultending+"_refine.txt"

    

    # postanalysis.writeedgetofile(outputfilepath1,postanalysis.removeunwantedgefromfile(resultweightfile1,importancethreshold=1.3))
    # postanalysis.writeedgetofile(outputfilepath2,postanalysis.removeunwantedgefromfile(resultweightfile2,importancethreshold=1.3))

    # refinedir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/testresult/"
    # outputfilepath= "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/testresult.txt"
    # postanalysis.findoutCommonedge(refinedir,outputfilepath)
            
            
    ########################################################################################
    # Post analysis
    #     
    ########################################################################################

    ########################################################################################
    # Whole process
    #     
    ########################################################################################
    # israndom = True
    # importancethreshold = 1.3
    # randomfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/new_result_pvalue_withppi_tyrobp.txt"
    # print("loading random dict")
    # randomnormaldic = postanalysis.loadrandomfile(randomfilepath)
    # inputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/"
    # outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/GSE115469_allcell/"
    # if not os.path.exists(outputdir):
    #     os.makedirs(outputdir)

    # print("loading graph")
    # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/tyrobp_related/sourcelink_tyrobp.txt'
    # sourcelist = kgmlparsar.loadjsontocontent(jsonfile)
    # iter = 2
    # subgroupnumber = 4

    
    # files= os.listdir(inputdir)
    # for file in files: 
    #     if not os.path.isdir(file): 
    #         refinefilepath = inputdir+file
    #         print("loading expression file "+file)
    #         allsamples = geneexpressionprocess.loadtotalfile(refinefilepath)  
    #         weight={}
    #         for i in range(iter):
    #             print("iteration "+str(i))
    #             # sperate
    #             subgroupsample = geneexpressionprocess.seperateall(subgroupnumber,allsamples,israndom)                
    #             geneindex= 0
    #             for targetgeneid,featuregeneids in sourcelist.items():                             
    #                 if targetgeneid in gene_id_symbol_map.keys():
    #                     targetgene = gene_id_symbol_map[targetgeneid]
    #                     geneindex+=1  
    #                     # if geneindex%100 ==0:
    #                     print(("gene id "+str(geneindex))+" "+targetgene +" : "+str(len(sourcelist)))
    #                     if not targetgene  in weight.keys():
    #                         weight[targetgene]={}
                        
    #                     featuregenelist = []
    #                     for nodeid in featuregeneids:
    #                         if nodeid in gene_id_symbol_map.keys():
    #                             featuregenelist.append(gene_id_symbol_map[nodeid])            
    #                     # for testgroup1
    #                     # index = 0
    #                     for samples in subgroupsample:
    #                         # index+=1
    #                         # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
    #                         weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #                         for featurename,weightscore in weightlist.items():
    #                             if not featurename in weight[targetgene].keys():
    #                                 weight[targetgene][featurename]=[weightscore]
    #                             else:
    #                                 weight[targetgene][featurename].append(weightscore)
            
    #         resultweightfilepath = outputdir+file                
    #         with open(resultweightfilepath,'w') as resultweightfile:
    #             for target,featureweights in weight.items():
    #                 for featurename,featureweight in featureweights.items():
    #                         pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                         resultweightfile.write(pline)
            
    #         postanalysis.reviseresultfile(resultweightfilepath,randomnormaldic)


            # postanalysis.writeedgetofile(resultweightfilepath,postanalysis.removeunwantedgefromfile(resultweightfilepath,importancethreshold))


    # commonedgefilepath = outputdir+"commonedge.txt"
    # postanalysis.findoutCommonedge(outputdir,commonedgefilepath)

    # os.system("shutdown -s -t  1")
    ########################################################################################
    # Whole process for specific genelist
    #     
    ########################################################################################
    # israndom = True
    # # importancethreshold = 1.3
    # # randomfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/new_result_pvalue_withppi_tyrobp.txt"
    # print("loading random dict")
    # # randomnormaldic = postanalysis.loadrandomfile(randomfilepath)
    # inputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/combinedranksum/0.1_small/"
    # outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/combinedranksum_small_high_6-13/"
    # if not os.path.exists(outputdir):
    #     os.makedirs(outputdir)
    # gene_symbol_id_map = kgmlparsar.getnametoid()
    # print("loading graph")
    # # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/tyrobp_related/sourcelink_tyrobp.txt'

    # checkgenelist=['AIF1','XBP1','GPX1','S100A9','S100A8','CTSS','LST1','SAT1','APOE']
    # # checkgenelist=['S100A8']
    # # checkgenelist=['IL6','TNF']
    # genelistset = kgmlparsar.getgeneidfromlist(checkgenelist,gene_symbol_id_map)

    # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allsourcelink_new_allgene_test_highconfidence.txt'
    # sourcelist = kgmlparsar.loadjsontocontent(jsonfile)
    
    # sourcelist = kgmlparsar.getsubsetgene(sourcelist, genelistset,model = 3)
    # print(len(sourcelist))
    # iter = 1
    # subgroupnumber = 4

    
    # files= os.listdir(inputdir)
    # for file in files: 
    #     if not os.path.isdir(file): 
    #         refinefilepath = inputdir+file
    #         print("loading expression file "+file)
    #         allsamples = geneexpressionprocess.loadtotalfile(refinefilepath)  
    #         weight={}
    #         for i in range(iter):
    #             print("iteration "+str(i))
    #             # sperate
    #             subgroupsample = geneexpressionprocess.seperateall(subgroupnumber,allsamples,israndom)                
    #             geneindex= 0
    #             for targetgeneid,featuregeneids in sourcelist.items():                             
    #                 if targetgeneid in gene_id_symbol_map.keys():
    #                     targetgene = gene_id_symbol_map[targetgeneid]
    #                     geneindex+=1  
    #                     # if geneindex%100 ==0:
    #                     print(("gene id "+str(geneindex))+" "+targetgene +" : "+str(len(sourcelist)))
    #                     if not targetgene  in weight.keys():
    #                         weight[targetgene]={}
                        
    #                     featuregenelist = []
    #                     for nodeid in featuregeneids:
    #                         if nodeid in gene_id_symbol_map.keys():
    #                             featuregenelist.append(gene_id_symbol_map[nodeid])            
    #                     # for testgroup1
    #                     # index = 0
    #                     for samples in subgroupsample:
    #                         # index+=1
    #                         # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
    #                         weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #                         for featurename,weightscore in weightlist.items():
    #                             if not featurename in weight[targetgene].keys():
    #                                 weight[targetgene][featurename]=[weightscore]
    #                             else:
    #                                 weight[targetgene][featurename].append(weightscore)
            
    #         resultweightfilepath = outputdir+file                
    #         with open(resultweightfilepath,'w') as resultweightfile:
    #             for target,featureweights in weight.items():
    #                 for featurename,featureweight in featureweights.items():
    #                         pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                         resultweightfile.write(pline)
            
            # postanalysis.reviseresultfile(resultweightfilepath,randomnormaldic)


            # postanalysis.writeedgetofile(resultweightfilepath,postanalysis.removeunwantedgefromfile(resultweightfilepath,importancethreshold))


    # commonedgefilepath = outputdir+"commonedge.txt"
    # postanalysis.findoutCommonedge(outputdir,commonedgefilepath)

    # os.system("shutdown -s -t  1")

#########################################################################
# compare all file
#
#########################################################################

    # resultfilelist1 = ["C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_cDC_2wks/commonedge_MCDD_Liver_c57bl_cDC_2wks.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_cDC_2wks/GSM5203271_AB2641.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_cDC_2wks/GSM5203272_AB2642.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_cDC_2wks/GSM5203275_AB795.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_cDC_2wks/GSM5203276_AB796.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_preDC_2wks/commonedge_MCDD_Liver_c57bl_preDC_2wks.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_preDC_2wks/GSM5203279_AB799.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_preDC_2wks/GSM5203280_AB800.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_cDC_2wks/commonedge_ND_Liver_c57bl_cDC_2wks.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_cDC_2wks/GSM5203269_AB2639.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_cDC_2wks/GSM5203270_AB2640.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_cDC_2wks/GSM5203273_AB793.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_cDC_2wks/GSM5203274_AB794.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_preDC_2wks/commonedge_ND_Liver_c57bl_preDC_2wks.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_preDC_2wks/GSM5203277_AB797.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_preDC_2wks/GSM5203278_AB798.txt",
    # ]
    # resultfilelist = [
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_cDC_2wks/commonedge_MCDD_Liver_c57bl_cDC_2wks.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/MCDD_Liver_c57bl_preDC_2wks/commonedge_MCDD_Liver_c57bl_preDC_2wks.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_cDC_2wks/commonedge_ND_Liver_c57bl_cDC_2wks.txt",
    # "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/ND_Liver_c57bl_preDC_2wks/commonedge_ND_Liver_c57bl_preDC_2wks.txt",
    # ]
    
    # outputfilepath = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/featureselectionsummarycommon.txt'
    # postanalysis.mappwithalledge(resultfilelist,outputfilepath)



    ########################################################################################
    # Replace the target value with mean value TEST
    # input control and test sample
    # not working 
    #     
    ########################################################################################
    # israndom = True
    # selectedsamplesize = 50
    # threshold=0.1
    # # importancethreshold = 1.3
    # # randomfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/new_result_pvalue_withppi_tyrobp.txt"
    # # print("loading random dict")
    # # randomnormaldic = postanalysis.loadrandomfile(randomfilepath)

    # controllabel =0
    # testlabel = 1
    # notshowinglabel = 2



    # # inputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/combinedranksum/0.1_small/"

    # controlsamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/Kupffer cell.txt"
    # testsamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/inflammatory macrophage.txt"


    # controlrefinesamples ="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/Kupffer cell_"+str(threshold)+".txt"
    # # "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl_refine/GSM5203286_control_DC_combined_"+str(threshold)+".txt"
    # testrefinesamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/inflammatory macrophage_"+str(threshold)+".txt"


    # preprocessdata.findthedifferentgene(controlsamples,testsamples,controlrefinesamples,testrefinesamples,threshold)


    # outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/GSE115469_refine/"
    # outputfilename="GSE115469_inflammatory_NORMAL_macrophage.txt"
    # if not os.path.exists(outputdir):
    #     os.makedirs(outputdir)
    # gene_symbol_id_map = kgmlparsar.getnametoid()
    # print("loading graph")
    # # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/tyrobp_related/sourcelink_tyrobp.txt'

    # checkgenelist=['AIF1','TYROBP','C1QA','C1QB','FTL','GPX1','LYZ','MARCO','S100A9','S100A8']
    # # checkgenelist=['S100A8']
    # # checkgenelist=['IL6','TNF']
    # genelistset = kgmlparsar.getgeneidfromlist(checkgenelist,gene_symbol_id_map)

    # # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allsourcelink_new_allgene_test.txt'
    # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allsourcelink_new_allgene_test_highconfidence.txt'
    # sourcelist = kgmlparsar.loadjsontocontent(jsonfile)
    
    # sourcelist = kgmlparsar.getsubsetgene(sourcelist, genelistset,model = 3)
    # print(len(sourcelist))
    # iter = 1
    # # subgroupnumber = 4

    
    # # files= os.listdir(inputdir)
    # # for file in files: 
    #     # if not os.path.isdir(file): 
    #         # refinefilepath = inputdir+file



    # # load samplefile        
    # print("loading "+ controlrefinesamples)
    # controlsamples = geneexpressionprocess.loadtotalfile(controlrefinesamples)

    # print("loading "+ testrefinesamples)
    # testsamples = geneexpressionprocess.loadtotalfile(testrefinesamples)  


    # weight={}
    # for i in range(iter):
    #     print("iteration "+str(i))
    #     # sperate
    #     controlsubgroupsample = geneexpressionprocess.randomselectedsamples(selectedsamplesize,controlsamples)
    #     testsubgroupsample = geneexpressionprocess.randomselectedsamples(selectedsamplesize,testsamples)

    #     # subgroupsample = geneexpressionprocess.seperateall(subgroupnumber,allsamples,israndom)                
    #     geneindex= 0

    #     for targetgeneid,featuregeneids in sourcelist.items():                             
    #         if targetgeneid in gene_id_symbol_map.keys():
    #             targetgene = gene_id_symbol_map[targetgeneid]
    #             geneindex+=1  
    #             # if geneindex%100 ==0:
    #             print(("gene id "+str(geneindex))+" "+targetgene +" : "+str(len(sourcelist)))
    #             if not targetgene  in weight.keys():
    #                 weight[targetgene]={}
                
    #             featuregenelist = []
    #             for nodeid in featuregeneids:
    #                 if nodeid in gene_id_symbol_map.keys():
    #                     featuregenelist.append(gene_id_symbol_map[nodeid])            
    #             # for testgroup1
    #             # index = 0

    #             # combine test and control as a samplegroupgroup
    #             totalsubgroup=[]
    #             targetspecificcontrolsubgroupsample=copy.deepcopy(controlsubgroupsample)
    #             targetspecifictestsubgroupsample=copy.deepcopy(testsubgroupsample)

    #             for sample in targetspecificcontrolsubgroupsample:
    #                 if targetgene in sample.keys():
    #                     sample[targetgene] = controllabel
    #                 else:
    #                     sample[targetgene] = notshowinglabel
    #                 totalsubgroup.append(sample)
                

    #             for sample in targetspecifictestsubgroupsample:
    #                 if targetgene in sample.keys():
    #                     sample[targetgene] = testlabel
    #                 else:
    #                     sample[targetgene] = notshowinglabel
    #                 totalsubgroup.append(sample)


    #             # for samples in subgroupsample:
    #                 # index+=1
    #                 # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
    #             weightlist = featureselectiongiveweight.processRelif(totalsubgroup,targetgene,featuregenelist)

    #             for featurename,weightscore in weightlist.items():
    #                 if not featurename in weight[targetgene].keys():
    #                     weight[targetgene][featurename]=[weightscore]
    #                 else:
    #                     weight[targetgene][featurename].append(weightscore)
            
    #         resultweightfilepath = outputdir+outputfilename                
    #         with open(resultweightfilepath,'w') as resultweightfile:
    #             for target,featureweights in weight.items():
    #                 for featurename,featureweight in featureweights.items():
    #                         pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                         resultweightfile.write(pline)
            
    #         # postanalysis.reviseresultfile(resultweightfilepath,randomnormaldic)


    #         # postanalysis.writeedgetofile(resultweightfilepath,postanalysis.removeunwantedgefromfile(resultweightfilepath,importancethreshold))


    # # commonedgefilepath = outputdir+"commonedge.txt"
    # # postanalysis.findoutCommonedge(outputdir,commonedgefilepath)

    # # os.system("shutdown -s -t  1")



# #########################################################################################################
# #
# try with feature selection on TRRUST
#
###########################################################################################################

    israndom = True
    selectedsamplesize = 200
    threshold=0.1
    # importancethreshold = 1.3
    # randomfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/new_result_pvalue_withppi_tyrobp.txt"
    # print("loading random dict")
    # randomnormaldic = postanalysis.loadrandomfile(randomfilepath)

    controllabel =0
    testlabel = 1
    notshowinglabel = 2



    # inputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/combinedranksum/0.1_small/"

    # testdatasample1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_1.txt"
    # testlabel1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_label1.txt"

    # testdatasample2 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_2.txt"
    # testlabel2 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_label2.txt"

    # testdatasample3 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_3.txt"
    # testlabel3 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_label3.txt"

    # testdatasample4 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_4.txt"
    # testlabel4 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_label4.txt"

    # testdatasample5 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_5.txt"
    # testlabel5 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE138852/split/AD_label5.txt"



    testdatasample1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control1.txt"
    testlabel1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control_label1.txt"

    testdatasample2 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control2.txt"
    testlabel2 =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control_label2.txt"


    testdatasample3 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control3.txt"
    testlabel3 =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control_label3.txt"


    testdatasample4 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control4.txt"
    testlabel4 =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control_label4.txt"


    testdatasample5 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control5.txt"
    testlabel5 =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/kmeanstransed_724/kmeans_data_control_label5.txt"







    # controlrefinesamples ="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/Kupffer cell_"+str(threshold)+".txt"
    # # "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl_refine/GSM5203286_control_DC_combined_"+str(threshold)+".txt"
    # testrefinesamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/inflammatory macrophage_"+str(threshold)+".txt"


    # preprocessdata.findthedifferentgene(controlsamples,testsamples,controlrefinesamples,testrefinesamples,threshold)

    outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-17/singlecellnormal/"

    # outputfilename="GSE115469_inflammatory_NORMAL_macrophage.txt"
    resultweightfile1 = outputdir+"random1.txt"
    resultweightfile2 = outputdir+"random2.txt"
    resultweightfile3 = outputdir+"random3.txt"
    resultweightfile4 = outputdir+"random4.txt"
    resultweightfile5 = outputdir+"random5.txt"
    # resultweightfile6 = outputdir+"kmeans_GSE169445_NASH_DC3.txt"
    # resultweightfile7 = outputdir+"kmeans_GSE169445_control_DC4.txt"
    # resultweightfile8 = outputdir+"kmeans_GSE169445_NASH_DC4.txt"
    # resultweightfile9 = outputdir+"kmeans_GSE169445_control_DC5.txt"
    # resultweightfile10 = outputdir+"kmeans_GSE169445_NASH_DC5.txt"
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    gene_symbol_id_map = kgmlparsar.getnametoid()
    print("loading graph")

    #########################################################################################################################################################################################################

    TFcol= 1
    target =2
    tfmapfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/TRRUST/trrust_rawdata.human.tsv"
    graphdict ={}

    with open(tfmapfilepath,"r") as graphfile:
        for line in graphfile.readlines():
            linedata= line.strip().split("\t")
            tfgenename = linedata[TFcol-1]
            targetgenename = linedata[target-1]            

            if not targetgenename in graphdict.keys():
                graphdict[targetgenename]=[]

            if not tfgenename in graphdict[targetgenename] and not tfgenename==targetgenename:
                graphdict[targetgenename].append(tfgenename)

    
		
    #########################################################################################################################################################################################################
    # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/tyrobp_related/sourcelink_tyrobp.txt'

    # checkgenelist=['AIF1','TYROBP','C1QA','C1QB','FTL','GPX1','LYZ','MARCO','S100A9','S100A8']
    # # checkgenelist=['S100A8']
    # # checkgenelist=['IL6','TNF']
    # genelistset = kgmlparsar.getgeneidfromlist(checkgenelist,gene_symbol_id_map)

    # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allsourcelink_new_allgene_test.txt'
    # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allsourcelink_new_allgene_test_highconfidence.txt'

    sourcelist = graphdict
    
    # sourcelist = kgmlparsar.getsubsetgene(sourcelist, genelistset,model = 3)

    print(len(sourcelist))
    iter = 1
    # subgroupnumber1 = 3
    # subgroupnumber2 = 5

    
    # files= os.listdir(inputdir)
    # for file in files: 
        # if not os.path.isdir(file): 
            # refinefilepath = inputdir+file



    # load samplefile        
    print("loading 1")
    dataset1 = geneexpressionprocess.loadtotalfile(testdatasample1)
    labelset1 = geneexpressionprocess.loadtotalfile(testlabel1)

    print("loading 2")
    dataset2 = geneexpressionprocess.loadtotalfile(testdatasample2)
    labelset2 = geneexpressionprocess.loadtotalfile(testlabel2)

    print("loading 3")
    dataset3 = geneexpressionprocess.loadtotalfile(testdatasample3)
    labelset3 = geneexpressionprocess.loadtotalfile(testlabel3)

    print("loading 4")
    dataset4 = geneexpressionprocess.loadtotalfile(testdatasample4)
    labelset4 = geneexpressionprocess.loadtotalfile(testlabel4)

    print("loading 5")
    dataset5 = geneexpressionprocess.loadtotalfile(testdatasample5)
    labelset5 = geneexpressionprocess.loadtotalfile(testlabel5)

   
    


    weight1={}
    weight2={}
    weight3={}
    weight4={}
    weight5={}
    

    for i in range(iter):
        print("iteration "+str(i))
        # sperate
        # subgroupsample = geneexpressionprocess.
        # controlsubgroupsamples1 = geneexpressionprocess.seperateall(subgroupnumber1,controlsamples1,israndom)  
        # testsubgroupsamples1 = geneexpressionprocess.seperateall(subgroupnumber2,testsamples1,israndom)  

        # controlsubgroupsamples2 = geneexpressionprocess.seperateall(subgroupnumber1,controlsamples2,israndom)  
        # testsubgroupsamples2 = geneexpressionprocess.seperateall(subgroupnumber2,testsamples2,israndom)  

        # controlsubgroupsamples3 = geneexpressionprocess.seperateall(subgroupnumber1,controlsamples3,israndom)  
        # testsubgroupsamples3 = geneexpressionprocess.seperateall(subgroupnumber2,testsamples3,israndom)  


        # subgroupsample = geneexpressionprocess.seperateall(subgroupnumber,allsamples,israndom)                
        geneindex= 0

        for targetgeneid,featuregeneids in sourcelist.items():                             
            # if targetgeneid in gene_id_symbol_map.keys():
            targetgene = targetgeneid
        # targetgene = 'NFKB1'
            geneindex+=1  
            # if geneindex%100 ==0:
            print(("gene id "+str(geneindex))+" "+targetgene +" : "+str(len(sourcelist)))
            if not targetgene  in weight1.keys():
                weight1[targetgene]={}
            
            if not targetgene  in weight2.keys():
                weight2[targetgene]={}

            if not targetgene  in weight3.keys():
                weight3[targetgene]={}
            
            if not targetgene  in weight4.keys():
                weight4[targetgene]={}
            
            if not targetgene  in weight5.keys():
                weight5[targetgene]={}
            
            # if not targetgene  in weight6.keys():
            #     weight6[targetgene]={}
            
            # if not targetgene  in weight7.keys():
            #     weight7[targetgene]={}
            
            # if not targetgene  in weight8.keys():
            #     weight8[targetgene]={}
            
            # if not targetgene  in weight9.keys():
            #     weight9[targetgene]={}
            
            # if not targetgene  in weight10.keys():
            #     weight10[targetgene]={}
            
            featuregenelist = sourcelist[targetgene]
        # featuregenelist = featuregeneids
        # for nodeid in featuregeneids:
        #     # if nodeid in gene_id_symbol_map.keys():
        #     featuregenelist.append(nodeid)            
        # for testgroup1
        # index = 0

        # combine test and control as a samplegroupgroup
        # totalsubgroup=[]
        # targetspecificcontrolsubgroupsample=copy.deepcopy(controlsubgroupsample)
        # targetspecifictestsubgroupsample=copy.deepcopy(testsubgroupsample)



        # for testgroup1
        # index = 0
            # for samples in controlsubgroupsamples1:
                # index+=1
                # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))




            weightlist = featureselectiongiveweight.processRelifwithlabel(dataset1,labelset1,targetgene,featuregenelist)
            for featurename,weight in weightlist.items():
                if not featurename in weight1[targetgene].keys():
                    weight1[targetgene][featurename]=[weight]
                else:
                    weight1[targetgene][featurename].append(weight)
            # index = 0
            # for samples in testsubgroupsamples1:
                # index+=1
                # print("subgroupnumber2 "+str(index)+":"+str(subgroupnumber))

            weightlist = featureselectiongiveweight.processRelifwithlabel(dataset2,labelset2,targetgene,featuregenelist)
            for featurename,weight in weightlist.items():
                if not featurename in weight2[targetgene].keys():
                    weight2[targetgene][featurename]=[weight]
                else:
                    weight2[targetgene][featurename].append(weight)
            
            # for samples in controlsubgroupsamples2:
                # index+=1
                # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
            weightlist = featureselectiongiveweight.processRelifwithlabel(dataset3,labelset3,targetgene,featuregenelist)
            for featurename,weight in weightlist.items():
                if not featurename in weight3[targetgene].keys():
                    weight3[targetgene][featurename]=[weight]
                else:
                    weight3[targetgene][featurename].append(weight)
            # index = 0
        # for samples in testsubgroupsamples2:
                # index+=1
                # print("subgroupnumber2 "+str(index)+":"+str(subgroupnumber))
            weightlist = featureselectiongiveweight.processRelifwithlabel(dataset4,labelset4,targetgene,featuregenelist)
            for featurename,weight in weightlist.items():
                if not featurename in weight4[targetgene].keys():
                    weight4[targetgene][featurename]=[weight]
                else:
                    weight4[targetgene][featurename].append(weight)

        # for samples in controlsubgroupsamples3:
                # index+=1
                # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
            weightlist = featureselectiongiveweight.processRelifwithlabel(dataset5,labelset5,targetgene,featuregenelist)
            for featurename,weight in weightlist.items():
                if not featurename in weight5[targetgene].keys():
                    weight5[targetgene][featurename]=[weight]
                else:
                    weight5[targetgene][featurename].append(weight)
            # index = 0
            # for samples in testsubgroupsamples3:
                # index+=1
                # print("subgroupnumber2 "+str(index)+":"+str(subgroupnumber))
            # weightlist = featureselectiongiveweight.processRelif(testsamples3,targetgene,featuregenelist)
            # for featurename,weight in weightlist.items():
            #     if not featurename in weight6[targetgene].keys():
            #         weight6[targetgene][featurename]=[weight]
            #     else:
            #         weight6[targetgene][featurename].append(weight)    



            #   # for samples in controlsubgroupsamples3:
            #     # index+=1
            #     # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
            # weightlist = featureselectiongiveweight.processRelif(controlsamples4,targetgene,featuregenelist)
            # for featurename,weight in weightlist.items():
            #     if not featurename in weight7[targetgene].keys():
            #         weight7[targetgene][featurename]=[weight]
            #     else:
            #         weight7[targetgene][featurename].append(weight)
            # # index = 0
            # # for samples in testsubgroupsamples3:
            #     # index+=1
            #     # print("subgroupnumber2 "+str(index)+":"+str(subgroupnumber))
            # weightlist = featureselectiongiveweight.processRelif(testsamples4,targetgene,featuregenelist)
            # for featurename,weight in weightlist.items():
            #     if not featurename in weight8[targetgene].keys():
            #         weight8[targetgene][featurename]=[weight]
            #     else:
            #         weight8[targetgene][featurename].append(weight)  


            # # for samples in controlsubgroupsamples3:
            #     # index+=1
            #     # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
            # weightlist = featureselectiongiveweight.processRelif(controlsamples5,targetgene,featuregenelist)
            # for featurename,weight in weightlist.items():
            #     if not featurename in weight9[targetgene].keys():
            #         weight9[targetgene][featurename]=[weight]
            #     else:
            #         weight9[targetgene][featurename].append(weight)
            # # index = 0
            # # for samples in testsubgroupsamples3:
            #     # index+=1
            #     # print("subgroupnumber2 "+str(index)+":"+str(subgroupnumber))
            # weightlist = featureselectiongiveweight.processRelif(testsamples5,targetgene,featuregenelist)
            # for featurename,weight in weightlist.items():
            #     if not featurename in weight10[targetgene].keys():
            #         weight10[targetgene][featurename]=[weight]
            #     else:
            #         weight10[targetgene][featurename].append(weight)  




          
            
    # after process
    # need to figure out later, right now just calculate mean of each edge
    with open(resultweightfile1,'w') as resultweightfile:
        for target,featureweights in weight1.items():
            for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
                # if statistics.mean(featureweight) >= importancethreshold:
                    pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
                    resultweightfile.write(pline)
    
    with open(resultweightfile2,'w') as resultweightfile:
        for target,featureweights in weight2.items():
            for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
                # if statistics.mean(featureweight) >= importancethreshold:
                    pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
                    resultweightfile.write(pline)
    
    with open(resultweightfile3,'w') as resultweightfile:
        for target,featureweights in weight3.items():
            for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
                # if statistics.mean(featureweight) >= importancethreshold:
                    pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
                    resultweightfile.write(pline)
    
    with open(resultweightfile4,'w') as resultweightfile:
        for target,featureweights in weight4.items():
            for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
                # if statistics.mean(featureweight) >= importancethreshold:
                    pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
                    resultweightfile.write(pline)
    
    with open(resultweightfile5,'w') as resultweightfile:
        for target,featureweights in weight5.items():
            for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
                # if statistics.mean(featureweight) >= importancethreshold:
                    pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
                    resultweightfile.write(pline)
    
    # with open(resultweightfile6,'w') as resultweightfile:
    #     for target,featureweights in weight6.items():
    #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    #             # if statistics.mean(featureweight) >= importancethreshold:
    #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                 resultweightfile.write(pline)
    

    # with open(resultweightfile7,'w') as resultweightfile:
    #     for target,featureweights in weight7.items():
    #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    #             # if statistics.mean(featureweight) >= importancethreshold:
    #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                 resultweightfile.write(pline)
    
    # with open(resultweightfile8,'w') as resultweightfile:
    #     for target,featureweights in weight8.items():
    #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    #             # if statistics.mean(featureweight) >= importancethreshold:
    #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                 resultweightfile.write(pline)
    

    # with open(resultweightfile9,'w') as resultweightfile:
    #     for target,featureweights in weight9.items():
    #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    #             # if statistics.mean(featureweight) >= importancethreshold:
    #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                 resultweightfile.write(pline)

    # with open(resultweightfile10,'w') as resultweightfile:
    #     for target,featureweights in weight10.items():
    #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    #             # if statistics.mean(featureweight) >= importancethreshold:
    #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                 resultweightfile.write(pline)
            # postanalysis.reviseresultfile(resultweightfilepath,randomnormaldic)


            # postanalysis.writeedgetofile(resultweightfilepath,postanalysis.removeunwantedgefromfile(resultweightfilepath,importancethreshold))


    # commonedgefilepath = outputdir+"commonedge.txt"
    # postanalysis.findoutCommonedge(outputdir,commonedgefilepath)

    # os.system("shutdown -s -t  1")





# #########################################################################################################################################################################################################




    # israndom = True
    # selectedsamplesize = 200
    # threshold=0.1
    # # importancethreshold = 1.3
    # # randomfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Pyscript/Findedgeweight/new_result_pvalue_withppi_tyrobp.txt"
    # # print("loading random dict")
    # # randomnormaldic = postanalysis.loadrandomfile(randomfilepath)

    # controllabel =0
    # testlabel = 1
    # notshowinglabel = 2



    # # inputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/combinedranksum/0.1_small/"

    # controlsamples1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE169445_control_DC_combined1_small.txt"
    # # testsamples1 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE169445_NASH_DC_combined1.txt"

    # # controlsamples2 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE169445_control_DC_combined2.txt"
    # # testsamples2 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE169445_NASH_DC_combined2.txt"

    # # controlsamples3 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE169445_control_DC_combined3.txt"
    # # testsamples3 = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE169445_NASH_DC_combined3.txt"

    # # controlrefinesamples ="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/Kupffer cell_"+str(threshold)+".txt"
    # # # "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl_refine/GSM5203286_control_DC_combined_"+str(threshold)+".txt"
    # # testrefinesamples = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/inflammatory macrophage_"+str(threshold)+".txt"


    # # preprocessdata.findthedifferentgene(controlsamples,testsamples,controlrefinesamples,testrefinesamples,threshold)

    # outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_small/"

    # # outputfilename="GSE115469_inflammatory_NORMAL_macrophage.txt"
    # resultweightfile1 = outputdir+"GSE169445_control_DC1_small.txt"
    # # resultweightfile2 = outputdir+"GSE169445_NASH_DC1.txt"
    # # resultweightfile3 = outputdir+"GSE169445_control_DC2.txt"
    # # resultweightfile4 = outputdir+"GSE169445_NASH_DC2.txt"
    # # resultweightfile5 = outputdir+"GSE169445_control_DC3.txt"
    # # resultweightfile6 = outputdir+"GSE169445_NASH_DC3.txt"
    # if not os.path.exists(outputdir):
    #     os.makedirs(outputdir)
    # gene_symbol_id_map = kgmlparsar.getnametoid()
    # print("loading graph")

    # #########################################################################################################################################################################################################

    # TFcol= 1
    # target =2
    # tfmapfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/TRRUST/trrust_rawdata.human.tsv"
    # graphdict ={}

    # with open(tfmapfilepath,"r") as graphfile:
    #     for line in graphfile.readlines():
    #         linedata= line.strip().split("\t")
    #         tfgenename = linedata[TFcol-1]
    #         targetgenename = linedata[target-1]
            

    #         if not targetgenename in graphdict.keys():
    #             graphdict[targetgenename]=[]

    #         if not tfgenename in graphdict[targetgenename] and not tfgenename==targetgenename:
    #             graphdict[targetgenename].append(tfgenename)

    
		
    # #########################################################################################################################################################################################################
    # # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/tyrobp_related/sourcelink_tyrobp.txt'

    # # checkgenelist=['AIF1','TYROBP','C1QA','C1QB','FTL','GPX1','LYZ','MARCO','S100A9','S100A8']
    # # # checkgenelist=['S100A8']
    # # # checkgenelist=['IL6','TNF']
    # # genelistset = kgmlparsar.getgeneidfromlist(checkgenelist,gene_symbol_id_map)

    # # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allsourcelink_new_allgene_test.txt'
    # # jsonfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/allsourcelink_new_allgene_test_highconfidence.txt'

    # sourcelist = graphdict
    
    # # sourcelist = kgmlparsar.getsubsetgene(sourcelist, genelistset,model = 3)

    # print(len(sourcelist))
    # iter = 1
    # # subgroupnumber1 = 3
    # # subgroupnumber2 = 5

    
    # # files= os.listdir(inputdir)
    # # for file in files: 
    #     # if not os.path.isdir(file): 
    #         # refinefilepath = inputdir+file



    # # load samplefile        
    # print("loading "+ controlsamples1)
    # controlsamples1 = geneexpressionprocess.loadtotalfile(controlsamples1)

    # # print("loading "+ testsamples1)
    # # testsamples1 = geneexpressionprocess.loadtotalfile(testsamples1)  

    # # print("loading "+ controlsamples2)
    # # controlsamples2 = geneexpressionprocess.loadtotalfile(controlsamples2)

    # # print("loading "+ testsamples2)
    # # testsamples2 = geneexpressionprocess.loadtotalfile(testsamples2) 

    # # print("loading "+ controlsamples3)
    # # controlsamples3 = geneexpressionprocess.loadtotalfile(controlsamples3)

    # # print("loading "+ testsamples3)
    # # testsamples3 = geneexpressionprocess.loadtotalfile(testsamples3) 
    


    # weight1={}
    # # weight2={}
    # # weight3={}
    # # weight4={}
    # # weight5={}
    # # weight6={}

    # for i in range(iter):
    #     print("iteration "+str(i))
    #     # sperate
    #     # subgroupsample = geneexpressionprocess.
    #     # controlsubgroupsamples1 = geneexpressionprocess.seperateall(subgroupnumber1,controlsamples1,israndom)  
    #     # testsubgroupsamples1 = geneexpressionprocess.seperateall(subgroupnumber2,testsamples1,israndom)  

    #     # controlsubgroupsamples2 = geneexpressionprocess.seperateall(subgroupnumber1,controlsamples2,israndom)  
    #     # testsubgroupsamples2 = geneexpressionprocess.seperateall(subgroupnumber2,testsamples2,israndom)  

    #     # controlsubgroupsamples3 = geneexpressionprocess.seperateall(subgroupnumber1,controlsamples3,israndom)  
    #     # testsubgroupsamples3 = geneexpressionprocess.seperateall(subgroupnumber2,testsamples3,israndom)  


    #     # subgroupsample = geneexpressionprocess.seperateall(subgroupnumber,allsamples,israndom)                
    #     geneindex= 0

    #     for targetgeneid,featuregeneids in sourcelist.items():                             
    #         # if targetgeneid in gene_id_symbol_map.keys():
    #         targetgene = targetgeneid
    #     # targetgene = 'NFKB1'
    #         geneindex+=1  
    #         # if geneindex%100 ==0:
    #         print(("gene id "+str(geneindex))+" "+targetgene +" : "+str(len(sourcelist)))
    #         if not targetgene  in weight1.keys():
    #             weight1[targetgene]={}
            
    #         # if not targetgene  in weight2.keys():
    #         #     weight2[targetgene]={}

    #         # if not targetgene  in weight3.keys():
    #         #     weight3[targetgene]={}
            
    #         # if not targetgene  in weight4.keys():
    #         #     weight4[targetgene]={}
            
    #         # if not targetgene  in weight5.keys():
    #         #     weight5[targetgene]={}
            
    #         # if not targetgene  in weight6.keys():
    #         #     weight6[targetgene]={}
            
    #         featuregenelist = sourcelist[targetgene]
    #     # featuregenelist = featuregeneids
    #     # for nodeid in featuregeneids:
    #     #     # if nodeid in gene_id_symbol_map.keys():
    #     #     featuregenelist.append(nodeid)            
    #     # for testgroup1
    #     # index = 0

    #     # combine test and control as a samplegroupgroup
    #     # totalsubgroup=[]
    #     # targetspecificcontrolsubgroupsample=copy.deepcopy(controlsubgroupsample)
    #     # targetspecifictestsubgroupsample=copy.deepcopy(testsubgroupsample)



    #     # for testgroup1
    #     # index = 0
    #         # for samples in controlsamples1:
    #         # index+=1
    #         # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
    #         weightlist = featureselectiongiveweight.processRelif(controlsamples1,targetgene,featuregenelist)
    #         for featurename,weight in weightlist.items():
    #             if not featurename in weight1[targetgene].keys():
    #                 weight1[targetgene][featurename]=[weight]
    #             else:
    #                 weight1[targetgene][featurename].append(weight)
    #         # index = 0
    #         # for samples in testsubgroupsamples1:
    #         #     # index+=1
    #         #     # print("subgroupnumber2 "+str(index)+":"+str(subgroupnumber))
    #         #     weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #         #     for featurename,weight in weightlist.items():
    #         #         if not featurename in weight2[targetgene].keys():
    #         #             weight2[targetgene][featurename]=[weight]
    #         #         else:
    #         #             weight2[targetgene][featurename].append(weight)
            
    #         # for samples in controlsubgroupsamples2:
    #         #     # index+=1
    #         #     # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
    #         #     weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #         #     for featurename,weight in weightlist.items():
    #         #         if not featurename in weight3[targetgene].keys():
    #         #             weight3[targetgene][featurename]=[weight]
    #         #         else:
    #         #             weight3[targetgene][featurename].append(weight)
    #         # # index = 0
    #         # for samples in testsubgroupsamples2:
    #         #     # index+=1
    #         #     # print("subgroupnumber2 "+str(index)+":"+str(subgroupnumber))
    #         #     weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #         #     for featurename,weight in weightlist.items():
    #         #         if not featurename in weight4[targetgene].keys():
    #         #             weight4[targetgene][featurename]=[weight]
    #         #         else:
    #         #             weight4[targetgene][featurename].append(weight)

    #         # for samples in controlsubgroupsamples3:
    #         #     # index+=1
    #         #     # print("subgroupnumber1 "+str(index)+":"+str(subgroupnumber))
    #         #     weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #         #     for featurename,weight in weightlist.items():
    #         #         if not featurename in weight5[targetgene].keys():
    #         #             weight5[targetgene][featurename]=[weight]
    #         #         else:
    #         #             weight5[targetgene][featurename].append(weight)
    #         # # index = 0
    #         # for samples in testsubgroupsamples3:
    #         #     # index+=1
    #         #     # print("subgroupnumber2 "+str(index)+":"+str(subgroupnumber))
    #         #     weightlist = featureselectiongiveweight.processRelif(samples,targetgene,featuregenelist)
    #         #     for featurename,weight in weightlist.items():
    #         #         if not featurename in weight6[targetgene].keys():
    #         #             weight6[targetgene][featurename]=[weight]
    #         #         else:
    #         #             weight6[targetgene][featurename].append(weight)    
    #         break



          
            
    # # after process
    # # need to figure out later, right now just calculate mean of each edge
    # with open(resultweightfile1,'w') as resultweightfile:
    #     for target,featureweights in weight1.items():
    #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    #             # if statistics.mean(featureweight) >= importancethreshold:
    #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    #                 resultweightfile.write(pline)
    
    # # with open(resultweightfile2,'w') as resultweightfile:
    # #     for target,featureweights in weight2.items():
    # #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    # #             # if statistics.mean(featureweight) >= importancethreshold:
    # #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    # #                 resultweightfile.write(pline)
    
    # # with open(resultweightfile3,'w') as resultweightfile:
    # #     for target,featureweights in weight3.items():
    # #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    # #             # if statistics.mean(featureweight) >= importancethreshold:
    # #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    # #                 resultweightfile.write(pline)
    
    # # with open(resultweightfile4,'w') as resultweightfile:
    # #     for target,featureweights in weight4.items():
    # #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    # #             # if statistics.mean(featureweight) >= importancethreshold:
    # #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    # #                 resultweightfile.write(pline)
    
    # # with open(resultweightfile5,'w') as resultweightfile:
    # #     for target,featureweights in weight5.items():
    # #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    # #             # if statistics.mean(featureweight) >= importancethreshold:
    # #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    # #                 resultweightfile.write(pline)
    
    # # with open(resultweightfile6,'w') as resultweightfile:
    # #     for target,featureweights in weight6.items():
    # #         for featurename,featureweight in dict(sorted(featureweights.items(), key=lambda item: item[1],reverse=True)).items():
    # #             # if statistics.mean(featureweight) >= importancethreshold:
    # #                 pline = target+"\t"+featurename+"\t"+str(statistics.mean(featureweight))+'\n'
    # #                 resultweightfile.write(pline)
            
    #         # postanalysis.reviseresultfile(resultweightfilepath,randomnormaldic)


    #         # postanalysis.writeedgetofile(resultweightfilepath,postanalysis.removeunwantedgefromfile(resultweightfilepath,importancethreshold))


    # # commonedgefilepath = outputdir+"commonedge.txt"
    # # postanalysis.findoutCommonedge(outputdir,commonedgefilepath)

    # os.system("shutdown -s -t  1")
