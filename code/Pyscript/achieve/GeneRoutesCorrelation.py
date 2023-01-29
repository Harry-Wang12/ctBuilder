

from scipy import stats
import itertools
from multiprocessing import Process



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

def calculateGeneRouteCorrelation(geneRoutespath,R1list,R2list,GeneCount,C,correlationType,importantGeneRouteRbar,importantGeneRoutePvalbar,rslist):
    with open(geneRoutespath,"w") as gengRoutefile:
        gengRoutefile.write("Gene\t"+C[0]+"\t"+C[1]+"\t"+str(rslist['cor'])+"\n")
        # calculateGeneCor
        # i = 0
        for genename,genevaluelist in GeneCount.items():
            # i+=1
            # process = i/len(GeneCount) * 100
            # print(geneRoutespath,":",str(round(process,4)),"%")
            rg1slist = calculatecorrelation_mult(correlationType,C[0], genename, R1list, genevaluelist)
            rg2slist = calculatecorrelation_mult(correlationType,C[1], genename, R2list, genevaluelist)
            # if abs(rg1slist['cor']) >= importantGeneRouteRbar and abs(rg2slist['cor']) >=importantGeneRouteRbar and rg1slist['pval']<=importantGeneRoutePvalbar and rg2slist['pval']<=importantGeneRoutePvalbar:
            # print(genename)
            pline = genename+"\t"+str(rg1slist['cor'])+"\t"+str(rg2slist['cor'])+"\n"
            gengRoutefile.write(pline)





if __name__ == "__main__":

    importantRouteRbar = 0.5
    importantGeneRouteRbar = 0.5

    importantRoutePvalbar =0.05
    importantGeneRoutePvalbar =0.05

    correlationType = "pearsonr"  #pearsonr ,  spearmanr, kendalltau

    routeScorepath =  "./HFD_pathways_route_score.txt"
    geneRatiopath = "./HFD_log2combine.txt"

    # geneRoutespath = "./HFD_geneRoutesCor.txt"

    RouteContent={}
    # read routeScore table
   
    with open (routeScorepath,"r") as routeScorefile:
        title = True
        for line in routeScorefile.readlines():
            if title:
                title =False
            else:
                linedata = line.strip().split("\t")
                routename = linedata[0].upper()
                del(linedata[0])
                valuecontent = [float(i) for i in linedata]
                if sum(valuecontent)>0 and "P2" not in routename:
                    RouteContent[routename] = valuecontent


    GeneCount={}
    # read routeScore table
    with open (geneRatiopath,"r") as geneLogfile:
        title = True
        for line in geneLogfile.readlines():
            if title:
                title =False
            else:
                linedata = line.strip().split("\t")
                genename = linedata[0].upper()
                del(linedata[0])
                valuecontent = [float(i) for i in linedata]
                GeneCount[genename] = valuecontent


    print("startCalculating.........")
    # Calcualte route score correlation
    Cs = itertools.combinations(RouteContent.keys(), 2)
    for C in  Cs:           
        xlist =RouteContent[C[0]]
        ylist =RouteContent[C[1]]            
        rslist = calculatecorrelation_mult(correlationType,C[0], C[1], xlist, ylist)
        # if abs(rslist['cor']) >= importantRouteRbar and rslist['pval']<=importantRoutePvalbar:
        R1list =RouteContent[C[0]]
        R2list =RouteContent[C[1]]
        print("Find high correlation routes: ", C[0], " and ",C[1] )
        geneRoutespath = "./generoute_all/"+C[0]+"_"+C[1]+".txt"
        p = Process(target=calculateGeneRouteCorrelation, args=(geneRoutespath,R1list,R2list,GeneCount,C,correlationType,importantGeneRouteRbar,importantGeneRoutePvalbar,rslist,))
        p.start()
            
