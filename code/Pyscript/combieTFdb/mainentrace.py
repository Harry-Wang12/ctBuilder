

import TFtargetnetworkrelated
import preprocessdata
import kgmlparsar
import Graphrelated
import os
import handlepathwayroute
import geneexpressionprocess

from itertools import permutations

from multiprocessing import Process



def dircombination(p,betweenroutedir,genelists,TFgraphdict,nodelength,TFgenelist):

    dirname = betweenroutedir + p[0]+'_'+p[1]+"/"
    allgenelistfile = betweenroutedir + p[0]+'_'+p[1]+"_combinelist.txt"
    print(dirname)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    allgenelist = []
    genelist1 = genelists[p[0]]
    genelist2 = genelists[p[1]]
    
    startandendcombination = handlepathwayroute.findallgenecombine(genelist1,genelist2)
    for combinationse in startandendcombination:
        # combinationse = startandendcombination[i+j]
        start = combinationse[0]
        end = combinationse[1]
        outputfilepath=dirname+start+'_'+end+".txt"
        if not start==end and start in TFgenelist and end in TFgenelist and not os.path.isfile(outputfilepath):
            allpaths = handlepathwayroute.findpathandwrite(TFgraphdict,start,end,nodelength,outputfilepath,genelist1,genelist2)
            for path in allpaths:
                for gene in path:
                    if not gene in allgenelist:
                        allgenelist.append(gene)
    with open(allgenelistfile,'w') as allgenelistoutput:
        for gene in allgenelist:
            allgenelistoutput.write(gene+'\n')

                # p = Process(target = handlepathwayroute.finepathandwrite,args=(TFgraphdict,start,end,nodelength,outputfilepath))
                # p.start()
                # l.append(p) 






if __name__ == '__main__':

    # filepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/TRRUST/trrust_rawdata.human.tsv"
    # start = 'GNA12'
    # end = 'ABL1'

    # # # filepath1="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/inflammatory macrophage.txt"
    # # # filepath2="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/Kupffer cell.txt"
    # # threshold=0.5
    # # # DEgeneoutputfilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE115469/allliver/GSE115469_inflammatory_control_"+str(threshold)+"_DE.txt"

    # # DEgeneoutputfilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE169445/NASHandControl/GSE5203285_NASH_control_"+str(threshold)+"_DE.txt"
    # # print("finding DEgene")
    # # # degenelist = preprocessdata.findthedifferentgene(filepath1,filepath2,DEgeneoutputfilepath,threshold)

    # # degenelist = TFtargetnetworkrelated.loadDElist(DEgeneoutputfilepath)


    # # print("load graph")
    # graphdict = TFtargetnetworkrelated.loadgraphfromfile(filepath)
    # # graphdict=TFtargetnetworkrelated.filterlink(degenelist,graphdict)
    # pathlength = 5
    # # print("find link from", start," to ",end)
    # allpath = TFtargetnetworkrelated.Check(graphdict,start,end,pathlength)
    # for p in allpath:
    #     print("All paths :%s"%p)

    # test line
    # filepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/TRRUST/trrust_rawdata.human.tsv"
    filepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-3/remainlink/testremain_0.8_0.1.txt"
    # filepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_6-27/controlremain.txt"

    TFgraphdict,TFgenelist = TFtargetnetworkrelated.loadgraphfromfile(filepath)


    gene_id_symbol_map=kgmlparsar.getgeneidsymbolmap()
    TF_genelist = kgmlparsar.getTFgenelist()
    receptor_gene_list = kgmlparsar.getreceptorlist()

    dirpath= 'C:/Users/whl19/Documents/Code/GenebetweenPathways/kgmlfiles_NASH/'

   

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
    



    # routefilepath = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/pathways/combined_route.txt'
    # outputdir ="C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/pathways/Keggpathwayroutegene/"
    # getgenefromroute(outputdir,routefilepath)
    nodelength = 15
    betweenroutedir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-3/betweenroute/"
    # betweenroutedir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_6-27/betweenroute_control/"
    
    numthreads = 16
    # allgenelistfilepath = getallgenelistfilename(outputdir)
    
    filecombination = list(permutations(kgmlfilelist, 2))
    
    # needpathwayname = 'hsa04620.xml'
    i=0
    while True:  
        l=[] 
        for j in range(numthreads):
            if i+j<len(filecombination):
                # if needpathwayname in filecombination[i+j]:
                # p = filecombination[i+j]
                p = Process(target = dircombination,args=(filecombination[i+j],betweenroutedir,genelists,TFgraphdict,nodelength,TFgenelist))

                # combinationse = startandendcombination[i+j]
                # start = combinationse[0]
                # end = combinationse[1]
                # outputfilepath=dirname+start+'_'+end+".txt"
                # if not start==end and start in TFgenelist and end in TFgenelist and not os.path.isfile(outputfilepath):
                # # finepathandwrite(graphdict,start,end,nodelength,outputfilepath)
                #     p = Process(target = handlepathwayroute.finepathandwrite,args=(TFgraphdict,start,end,nodelength,outputfilepath))
                p.start()
                l.append(p) 

        for p in l :
            p.join() 
        i+=numthreads
        if i>=len(filecombination):
            break



            dirname = betweenroutedir + p[0]+'_'+p[1]+"/"
            print(dirname)
            if not os.path.exists(dirname):
                os.makedirs(dirname)

                
            genelist1 = genelists[p[0]]
            genelist2 = genelists[p[1]]





    # for p in permutations(kgmlfilelist, 2):

    #     dirname = betweenroutedir + p[0]+'_'+p[1]+"/"

    #     print(dirname)
    #     if not os.path.exists(dirname):
    #         os.makedirs(dirname)

               
    #     genelist1 = genelists[p[0]]
    #     genelist2 = genelists[p[1]]

    #     startandendcombination = handlepathwayroute.findallgenecombine(genelist1,genelist2)
        
    #     # i=0
    #     # while True:  
    #     l=[] 
    #     #     for j in range(numthreads):
    #     #         if i+j<len(startandendcombination): 
    #     for combinationse in startandendcombination:
    #                 # combinationse = startandendcombination[i+j]
    #         start = combinationse[0]
    #         end = combinationse[1]
    #         outputfilepath=dirname+start+'_'+end+".txt"
    #         if not start==end and start in TFgenelist and end in TFgenelist and not os.path.isfile(outputfilepath):
    #         # finepathandwrite(graphdict,start,end,nodelength,outputfilepath)
    #             p = Process(target = handlepathwayroute.finepathandwrite,args=(TFgraphdict,start,end,nodelength,outputfilepath))
    #             p.start()
    #             l.append(p) 

    #     for p in l :
    #         p.join() 
    #         i+=numthreads
    #     if i>=len(startandendcombination):
    #         break



        # for combinationstartend in startandendcombination:
        #     subprocess




    # pathwaykgmlfile = dirpath+"hsa04010.xml"
    # kgmlcontent = kgmlparsar.readkgml(pathwaykgmlfile)
    # graphdict,genelist = kgmlparsar.generateGraph(kgmlcontent,gene_id_symbol_map)

    # p1paths = Graphrelated.findp1route(graphdict,genelist,TF_genelist,receptor_gene_list)
    # print(p1paths)



    # # init 
    # print("graph")
    # for path,dir_list,file_list in os.walk("C:/Users/whl19/Documents/Code/GenebetweenPathways/kgmlfiles/"):
    #     for pathwayname in file_list:
            
    #         print("graph "+pathwayname)
    #         betweenroutedir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_7-4/betweenroute_NASH/"
    #         kgmlfilepath =  'C:/Users/whl19/Documents/Code/GenebetweenPathways/kgmlfiles/'+pathwayname
    #         outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_7-4/upperpath/"
    #         realoutputdir = outputdir+pathwayname+'/'
    #         if not os.path.exists(realoutputdir):
    #             os.makedirs(realoutputdir)
    #         gene_id_symbol_map = kgmlparsar.getgeneidsymbolmap()
    #         ligand_gene_list = kgmlparsar.getligandlist()
    #         receptor_gene_list = kgmlparsar.getreceptorlist()
    #         upsteamends = Graphrelated.findtheendoftheuppersteam(pathwayname,betweenroutedir)
    #         kgmlcontent = kgmlparsar.readkgml(kgmlfilepath,gene_id_symbol_map)
    #         kgmlgraphdict,kgmlgenelist = kgmlparsar.generateGraph(kgmlcontent,gene_id_symbol_map)
    #         # print(kgmlgraphdict['AKT1'])
    #         # print(kgmlgraphdict['PIK3R2'])
    #         # print(kgmlgraphdict['RAC1'])
    #         # print(kgmlgraphdict['TLR1'])
    #         # print(kgmlgraphdict['TLR2'])
    #         for endnode in upsteamends:
    #             Graphrelated.findupsteam(endnode,realoutputdir,receptor_gene_list,kgmlgraphdict,kgmlgenelist,ligand_gene_list)


    # print("sample")
    # # scoring
    # log2sampledata = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/GSE44770/LOAD/GSE44770_mapped.txt"
    # print("loading "+ log2sampledata)
    # allsamples = geneexpressionprocess.loadtotalfile(log2sampledata)
    # allsplitsampls = geneexpressionprocess.seperateall(3,allsamples,True)
    
    # print("up")
    # # upstreamlist=[]
    
    # # g = os.walk(realoutputdir)
    # # for path,dir_list,file_list in g:  
    # #     for file_name in file_list:  
    # #         upstreamlist+=handlepathwayroute.getpathfromfile(realoutputdir+file_name)

    # # upstreampaths=[]


    # # for uppath in upstreamlist:
    # #     # upper
    # #         uppatterns = handlepathwayroute.getupperstreampattern(kgmlcontent['edgelist'],uppath)
    # #         for uppattern in uppatterns:
    # #             pathcombo = {
    # #                 'pathname': uppath,
    # #                 'pathpattern':uppattern
    # #             }
    # #             upstreampaths.append(pathcombo)
    # print("middle")
    # betweenroutedir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/pathways/betweenroute/"

    # midestreamlist=[]

    # g = os.walk(betweenroutedir)
    # for path,dir_list,file_list in g:  
    #     for dir_name in dir_list:
    #         if dir_name.startswith(pathwayname):
    #             upperstreamdir = betweenroutedir+dir_name+'/'
    #             q = os.walk(upperstreamdir) 
    #             for qpath,qdir_list,qfile_list in q:
    #                 for file_name in qfile_list:
    #                     qfilepath = upperstreamdir+file_name
    #                     if os.path.getsize(qfilepath)>1:                            
    #                         midestreamlist+=handlepathwayroute.getpathfromfile(qfilepath)

    
    # relationfilepath="C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/TRRUST/trrust_rawdata.human_refine.tsv"

    # midrelation = handlepathwayroute.getrealtion(relationfilepath)

    # midpaths = []

    # for midpath in midestreamlist:
    #     # upper
    #         midpatterns = handlepathwayroute.getmidperstreampattern(midrelation,midpath)
    #         for midpattern in midpatterns:
    #             pathcombo = {
    #                 'pathname': midpath,
    #                 'pathpattern':midpattern
    #             }
    #             midpaths.append(pathcombo)

    
    # results=[]
    
    # index=1
    # scoreoutputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/up_middle/"

    # if not os.path.exists(scoreoutputdir):
    #     os.makedirs(scoreoutputdir)

    # # outputfilename="GSE115469_inflammatory_NORMAL_macrophage.txt"
    # # resultweightfile1path = scoreoutputdir+"1.txt"
    # # resultweightfile2path = scoreoutputdir+"2.txt"
    # # resultweightfile3path = scoreoutputdir+"3.txt"
    # print("scoring")
    # for samples in allsplitsampls:
    #     result = {}
    #     samplenumber = 1
    #     for sample in samples:  
    #         print(samplenumber/(len(samples)))
    #         samplenumber+=1
    #         # index =1          
    #         # for uppath in upstreampaths:  
    #         #     # index+=1
    #         #     # print(index/(len(upstreampaths)+len(midpaths)))
    #         #     # print(uppath['pathname'])              
    #         #     pathname = ",".join(uppath['pathname'])+":"+",".join([str(x) for x in uppath['pathpattern']])
    #         #     if not pathname in result.keys():
    #         #         result[pathname] = [handlepathwayroute.scoring(uppath['pathname'],uppath['pathpattern'],sample)]
    #         #     else:
    #         #         result[pathname] .append(handlepathwayroute.scoring(uppath['pathname'],uppath['pathpattern'],sample))
    #         #     # break

    #         for midpath in midpaths: 
    #             # index+=1
    #             # print(index/(len(upstreampaths)+len(midpaths)))
    #             # print(midpath['pathname'])                     
    #             pathname = ",".join(midpath['pathname'])+":"+",".join([str(x) for x in midpath['pathpattern']])
    #             if not pathname in result.keys():
    #                 result[pathname] = [handlepathwayroute.scoring(midpath['pathname'],midpath['pathpattern'],sample)]
    #             else:
    #                 result[pathname] .append(handlepathwayroute.scoring(midpath['pathname'],midpath['pathpattern'],sample))
    #             # break
    #     with open(scoreoutputdir+str(index)+'_middle.txt','w') as scorefile:
            
    #         for pathname,scores in result.items():
    #             pline = pathname
    #             for score in scores:
    #                 pline+='\t'+str(score)
    #             scorefile.write(pline+'\n')
    #         index+=1

    # os.system("shutdown -s -t  1")    
