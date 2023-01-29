
import kgmlparsar
import os
from itertools import permutations

from shutil import copyfile






if __name__=="__main__":
    
  

    # dirpath= 'C:/Users/whl19/Documents/Code/GenebetweenPathways/kgmlfiles/'

    betweengenesetdir ="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-3/betweenroute/"

    outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-3/betweenroute_summary/"

    commonfilelist=[]
    copylist = []
    destlist = []


 
    # filecombinations=["hsa04024.xml_hsa04218.xml",
    #                 "hsa04062.xml_hsa04015.xml",
    #                 "hsa04062.xml_hsa04218.xml",
    #                 "hsa04062.xml_hsa04550.xml",
    #                 "hsa04064.xml_hsa04550.xml",
    #                 "hsa04066.xml_hsa04218.xml",
    #                 "hsa04066.xml_hsa04550.xml",
    #                 "hsa04071.xml_hsa04015.xml",
    #                 "hsa04071.xml_hsa04550.xml",
    #                 "hsa04110.xml_hsa04015.xml",
    #                 "hsa04110.xml_hsa04340.xml",
    #                 "hsa04110.xml_hsa04550.xml",
    #                 "hsa04151.xml_hsa04015.xml",
    #                 "hsa04151.xml_hsa04150.xml",
    #                 "hsa04151.xml_hsa04218.xml",
    #                 "hsa04151.xml_hsa04550.xml",
    #                 "hsa04151.xml_hsa04630.xml",
    #                 "hsa04151.xml_hsa04666.xml",
    #                 "hsa04151.xml_hsa04916.xml",
    #                 "hsa04151.xml_hsa04920.xml",
    #                 "hsa04217.xml_hsa04015.xml",
    #                 "hsa04217.xml_hsa04390.xml",
    #                 "hsa04217.xml_hsa04550.xml",
    #                 "hsa04217.xml_hsa04916.xml",
    #                 "hsa04218.xml_hsa04015.xml",
    #                 "hsa04218.xml_hsa04340.xml",
    #                 "hsa04218.xml_hsa04550.xml",
    #                 "hsa04218.xml_hsa04630.xml",
    #                 "hsa04218.xml_hsa04916.xml",
    #                 "hsa04218.xml_hsa04919.xml",
    #                 "hsa04330.xml_hsa04015.xml",
    #                 "hsa04330.xml_hsa04550.xml",
    #                 "hsa04350.xml_hsa04218.xml",
    #                 "hsa04350.xml_hsa04340.xml",
    #                 "hsa04350.xml_hsa04550.xml",
    #                 "hsa04390.xml_hsa04550.xml",
    #                 "hsa04510.xml_hsa04550.xml",
    #                 "hsa04550.xml_hsa04390.xml",
    #                 "hsa04630.xml_hsa04015.xml",
    #                 "hsa04630.xml_hsa04550.xml",
    #                 "hsa04660.xml_hsa04550.xml",
    #                 "hsa04668.xml_hsa04550.xml",
    #                 "hsa04916.xml_hsa04015.xml",
    #                 "hsa04917.xml_hsa04218.xml",
    #                 "hsa04917.xml_hsa04550.xml",
    #                 "hsa04919.xml_hsa04015.xml",
    #                 "hsa04919.xml_hsa04550.xml",
    #                 "hsa04920.xml_hsa04218.xml",
    #                 "hsa04920.xml_hsa04550.xml",
    #                 "hsa04922.xml_hsa04015.xml",
    #                 "hsa04922.xml_hsa04550.xml",
    #                 "hsa04926.xml_hsa04015.xml",
    #                 "hsa04926.xml_hsa04550.xml",
    #                 "hsa04927.xml_hsa04015.xml",
    #                 "hsa04927.xml_hsa04218.xml",
    #                 "hsa04927.xml_hsa04550.xml",
    #                 "hsa05120.xml_hsa04550.xml",
    #                 "hsa05226.xml_hsa04015.xml",
    #                 ]


    
    # outputfile="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-3/betweenroute_NASH_RELATED_10_allpath.txt"

    # for filecombination in filecombinations:

    g = os.walk(betweengenesetdir+"/")  

    index = 1
    
    for path,dirs,file_list in os.walk(betweengenesetdir+"/"):
        for dir_name in dirs:
            allpaths=[]
            for subpath,subdirs,subfile_list in os.walk(os.path.join(path,dir_name)):
                outputfile = outputdir+dir_name+'.txt'

                with open(outputfile,'w') as output:
                # for path,dirs,file_list in os.walk(betweengenesetdir):
                    for subfile in subfile_list:
                        if os.path.getsize(os.path.join(subpath,subfile))>0: 

                            with open(os.path.join(subpath,subfile),'r') as pathwayfile:
                                for line in pathwayfile.readlines():
                                    linepath = line.strip().split("\t")
                                    # del(linepath[0])
                                    # del(linepath[0])
                                    if not linepath in allpaths:
                                        pline = "path_"+str(index)+'\t'+"path_"+str(index)
                                        for gene in linepath:
                                            pline+="\t"+gene
                                        output.write(pline+'\n')
                                        print(index ,linepath)
                                        index +=1
                                        allpaths.append(linepath)




    # with open(outpuffile,'w') as output:


    #     index = 1
    #     for path in allpaths:
    #         pline = "path_"+str(index)+'\t'+"path_"+str(index)
    #         for gene in path:
    #             pline+="\t"+gene
    #         index +=1
            
    #         output.write(pline+'\n')

        

