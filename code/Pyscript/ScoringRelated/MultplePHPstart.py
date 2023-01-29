import subprocess
import sys
import os
from multiprocessing import Process




def callScoreingphpscript(inputtable,outputtable,pathwayfile):
    # print(pathwayfile+"           start")
    result = subprocess.run(
        ['C:/wamp64/bin/php/php5.6.40/php.exe', 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/score_pathways_cliforsinglepathway.php',inputtable,outputtable,pathwayfile],    # program and arguments
        stdout=subprocess.PIPE,  # capture stdout
        check=False               # raise exception if program fails
    )
    # print(pathwayfile+"           done") 



def multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir):

    if not os.path.exists(outputdir):
        os.makedirs(outputdir) 
    pathwaypathlist = []
    pathwaynamelist = []
    g = os.walk(revisedpathwaydir)
    for path,dir_list,file_list in g:  
        for file_name in file_list: 
            pathwaypathlist.append(os.path.join(path, file_name))
            pathwaynamelist.append(file_name)

    i=0
    while True:       
        l=[]
        for j in range(numthreads):
            if i+j<len(pathwaypathlist): 
                outputtable = outputdir+"summary_"+pathwaynamelist[i+j]
                # print('C:/wamp64/bin/php/php5.6.40/php.exe', 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/score_pathways_cliforsinglepathway.php',inputtable,outputtable,pathwaypathlist[i+j])
                p = Process(target = callScoreingphpscript,args=(inputtable,outputtable,pathwaypathlist[i+j]))
                p.start()
                l.append(p) 
        for p in l :
            p.join() 
        i+=numthreads
        if i>=len(pathwaynamelist):
            break


if __name__=="__main__":

    allresultdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory_test/"
    inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/inflammatory_macrophage/combined.txt"
    revisedpathwaydir=allresultdir+"Revisedpathway/"
    print(revisedpathwaydir)
    
    numthreads = 16
    outputdir=allresultdir+"revisedScoring/"
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)





    # allresultdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/"
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/non-inflammatory_macrophage/combined.txt"
    # revisedpathwaydir=allresultdir+"Revisedpathway/1/"
    # print(revisedpathwaydir)
    
    # numthreads = 16
    # outputdir=allresultdir+"revisedScoring/"
    # if not os.path.exists(outputdir):
    #     os.makedirs(outputdir)
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)


    # allresultdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/"
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/non-inflammatory_macrophage/combined.txt"
    # revisedpathwaydir=allresultdir+"Revisedpathway/2/"
    # print(revisedpathwaydir)
    
    # numthreads = 16
    # outputdir=allresultdir+"revisedScoring/"
    # if not os.path.exists(outputdir):
    #     os.makedirs(outputdir)
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)

    # allresultdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/"
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/non-inflammatory_macrophage/combined.txt"
    # revisedpathwaydir=allresultdir+"Revisedpathway/3/"
    # print(revisedpathwaydir)
    
    # numthreads = 16
    # outputdir=allresultdir+"revisedScoring/"
    # if not os.path.exists(outputdir):
    #     os.makedirs(outputdir)
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)

    # allresultdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/"
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/non-inflammatory_macrophage/combined.txt"
    # revisedpathwaydir=allresultdir+"Revisedpathway/4/"
    # print(revisedpathwaydir)
    
    # numthreads = 16
    # outputdir=allresultdir+"revisedScoring/"
    # if not os.path.exists(outputdir):
    #     os.makedirs(outputdir)
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)




    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-15-2021_GSE115469_inflamtory/Revisedpathway/2500/"
    # print(revisedpathwaydir)
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE155182/GSE155182_NASH_7m_liver_rescale.txt"
    # numthreads = 16
    # # outputdir=allresultdir+"revisedScoring/"
    # # if not os.path.exists(outputdir):
    # #     os.makedirs(outputdir)
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)

    # os.system("shutdown -s -t  1")

    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/Revisedpathway/2500low"
    # print(revisedpathwaydir)
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # numthreads = 16
    # outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/revisedScoring_test/"
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)
    
    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/Revisedpathway/5000low"
    # print(revisedpathwaydir)
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # numthreads = 16
    # outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/revisedScoring_test/"
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)

    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/Revisedpathway/7500low"
    # print(revisedpathwaydir)
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # numthreads = 16
    # outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/revisedScoring_test/"
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)

    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/Revisedpathway/10000low"
    # print(revisedpathwaydir)
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # numthreads = 16
    # outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/revisedScoring_test/"
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)

    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/Revisedpathway/15000low"
    # print(revisedpathwaydir)
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # numthreads = 16
    # outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/revisedScoring_test/"
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)



    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/Revisedpathway/rest"
    # print(revisedpathwaydir)
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # numthreads = 16
    # outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/revisedScoring_test/"
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)


    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/Revisedpathway/1000low"
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # numthreads = 16
    # outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/revisedScoring/"
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)

    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/Revisedpathway/10000low"
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # numthreads = 16
    # outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/revisedScoring/"
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)

    # revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/Revisedpathway/10000high"
    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Cancer_data/TCGA_SKCM/TCGA_SKCM_log2ratio.txt"
    # numthreads = 16
    # outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/2-28-2021/revisedScoring/"
    # multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)


    # for i in range(5):
    #     revisedpathwaydir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/Revisedpathway/"
    #     inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/singlecellsimulation/counts_"+str(i)+"_log2/combined.txt"
    #     numthreads = 16
    #     outputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/revisedScoring_counts_"+str(i)+"_log2/"
    #     multipleCallscoringscript(revisedpathwaydir,numthreads,inputtable,outputdir)


    # callScoreingphpscript(inputtable,outputtable,pathwayfile)





    # os.system("shutdown -s -t  1")

