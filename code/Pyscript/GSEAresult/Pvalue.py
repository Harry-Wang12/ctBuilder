import statistics
from scipy import stats
import numpy as np
import os


def loadfile(filepath):

    allsamples={}
    
    with open(filepath,'r') as allsamplefile:
        istitle=True
        for line in allsamplefile.readlines():
            linedata = line.strip().split("\t")
            if istitle:               
                istitle=False
            else:

                pathwayroutename = linedata[0].upper()+":"+linedata[1].upper()
                allsamples[pathwayroutename] = []
                # print(genename)
                del(linedata[0])
                del(linedata[0])
                del(linedata[-1])
                for genevalue in linedata:
                    allsamples[pathwayroutename].append(float(genevalue)) 
    
    return allsamples


def writefile(filepath,samples):

    with open(filepath,'w') as outputfile:
        for pathwayname,pvalues in samples.items():
            pline = pathwayname.replace(':','\t')
            for pvalue in pvalues:
                pline+='\t'+str(pvalue)
            outputfile.write(pline+'\n')



def calculatedPvalues(randomdataset, realdataset):
    returndataset = {}
    bw_method= 0.001
    upper_bound= 1.1
    for pathwayname, scores in realdataset.items():
        # print(pathwayname)
        returndataset[pathwayname]=[]

        # for pw_route_score in scores:
        #     number = 0
        #     for value in randomdataset[pathwayname]:
        #         if value >=pw_route_score:
        #             number+=1
            
        #     pval = number/len(randomdataset[pathwayname])
        #     if pval > 0.5:
        #         pval = 1 - pval





        kde = stats.gaussian_kde(randomdataset[pathwayname], bw_method=bw_method)

        for pw_route_score in scores:
            pval = kde.integrate_box_1d(pw_route_score, upper_bound)
            if pval > 0.5:
                pval = 1 - pval

            returndataset[pathwayname].append(pval)


    return returndataset






if __name__ == "__main__":

    for path,dir_list,file_list in os.walk("C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/Scoredroutes_withtreatpaper/"):
        for filename in file_list:
            print(filename)
            datasetfile = os.path.join(path,filename)

            referencefile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/Scoredroutes_random/"+filename
            # datasetfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_7-29/uniquepathconsistent_sign.txt"

            outputfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/Scoredroutespvalue/"+filename


            referencedataset = loadfile(referencefile)
            realdataset = loadfile(datasetfile)


            pvaluedataset = calculatedPvalues(referencedataset, realdataset)

            writefile(outputfile,pvaluedataset)




