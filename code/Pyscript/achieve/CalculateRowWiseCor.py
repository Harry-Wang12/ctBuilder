



from scipy import stats
import itertools
from multiprocessing import Pool
from multiprocessing import Manager
from multiprocessing import Process
import multiprocessing
import os

def calculatecorrelation_mult(correlationType,xname,yname,xlist,ylist):

    
    returndict={}
    returndict["xname"] = xname
    returndict["yname"] = yname
    if correlationType == "pearsonr":
        cor,pval = stats.pearsonr(xlist,ylist)
        returndict['cor'] = cor
        returndict['pval'] = pval
    elif correlationType == "spearmanr":
        cor,pval =stats.spearmanr(xlist,ylist)
        returndict['cor'] = cor
        returndict['pval'] = pval
    elif correlationType == "kendalltau":
        cor,pval = stats.kendalltau(xlist,ylist)
        returndict['cor'] = cor
        returndict['pval'] = pval
    
    return returndict




# calculate rowwise corelation 






if __name__ == "__main__":
    # manager = Manager()

    # pool = Pool(multiprocessing.cpu_count(),maxtasksperchild=1) 
    

    correlationType = "pearsonr"  #pearsonr ,  spearmanr, kendalltau

   
    

    # genelist = ["HDAC4"]
    genelist = "ALL"

    filepath="./Cancer_data/TCGA_COAD_pathways_route_score_score5_20201227_012736_7362_wo.txt"

    resultpath = "./Cancer_data/"+correlationType+"_TCGA_COAD_pathway_routes.txt"


    rownamelist = []
    tablecontent={}
    with open(filepath,"r") as tablefile:
        Header = True

        for line in tablefile.readlines():
            if Header:
                Header=False
            else:
                linedata = line.strip().split("\t")
                rowname = linedata[0].upper()
                del(linedata[0])
                
                valuecontent = [float(i) for i in linedata]
                # if sum(valuecontent)>0 and "P2" not in rowname:
                tablecontent[rowname]=valuecontent
                rownamelist.append(rowname)

    Cs = []
    if genelist =="ALL":
        genelist=rownamelist
        Cs = itertools.combinations(genelist, 2)
    elif len(genelist) ==1:
        Cs =[]
        for i in rownamelist:
            Cs.append([genelist[0],i])
    else:
        Cs = itertools.combinations(genelist, 2)   
        


    # write result:
    with  open (resultpath,"w") as resultfile:
        resultfile.write("X\tY\tcorrelation\tP-val\n")
        i = 0
        for C in  Cs:
            i+=1
            print(i)
            xlist =tablecontent[C[0]]
            ylist =tablecontent[C[1]]
            rslist = calculatecorrelation_mult(correlationType,C[0], C[1], xlist, ylist)
            if not str(rslist['cor']) == "nan":
                pline = rslist["xname"]+"\t"+ rslist["yname"]+"\t"+str(rslist['cor'])+"\t"+str(rslist['pval'])+"\n"
                resultfile.write(pline)


        # for rslist in return_list:
        #     pline = rslist["xname"]+"\t"+ rslist["yname"]+"\t"+str(rslist['cor'])+"\t"+str(rslist['pval'])+"\n"
        #     resultfile.write(pline)


    # os.system('shutdown -s')

    