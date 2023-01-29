import scipy.stats as st
import matplotlib.pyplot as plt

def getnoticep2(Summaryfilepath,outputresultfile,figname):

    AllCorrelationlist = []
    with open(Summaryfilepath,"r") as PathwayRoutefile:

        pathwaylists={}
        title = True
        for line in PathwayRoutefile.readlines():
            if title:
                title=False
            else:
                linedata = line.replace("\"","").strip().split("\t")
                pinformation = linedata[0].strip().split("]")[1].split("~")
                pathwayname = pinformation[0]
                routepart = pinformation[1]
                source = pinformation[2]
                target = pinformation[3]

                if not pathwayname in pathwaylists.keys():
                    pathwaylists[pathwayname]={
                        "p1":[],
                        "p2":[]
                    }
                valueintem = {
                    'source':source,
                    'target':target,
                    'valuelist':[],
                }
                for i in range(2,len(linedata)):
                    valueintem['valuelist'].append(float(linedata[i]))

                pathwaylists[pathwayname][routepart].append(valueintem)
    
    with open(outputresultfile,"w") as outputfile:
        for pathwayname,routes in pathwaylists.items():
            for p2route in routes["p2"]:
                title =pathwayname+"\tp2Source: "+p2route["source"]+",p2Target: "+p2route["target"]
                printline=""
                for p1route in routes["p1"]:
                    if p1route["target"] == p2route["source"]:
                        # find original
                        p1title = "p1Source: "+p1route["source"]+",p1Target: "+p1route["target"]
                        r,p =st.pearsonr(p2route["valuelist"], p1route["valuelist"]) 
                        AllCorrelationlist.append(r)
                        printline+=title+"\t"+p1title+"\t"+str(r)+"\t"+str(p)+"\n"
                if len(printline)==0:
                    printline+=title+"\n"
                outputfile.write(printline)

    kwargs = dict(histtype='stepfilled', alpha=0.3, bins=50)            
    plt.figure()
    plt.hist(AllCorrelationlist, **kwargs)
    plt.savefig(figname)
    plt.close()


if __name__=="__main__":
    pathwaysocrefile="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_inflamtory/RouteScore.txt"    
    P1p2fileoutput="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_inflamtory/OriginalCorrelation.txt"    
    figname = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_inflamtory/OriginalCorfig.jpg" 
    getnoticep2(pathwaysocrefile,P1p2fileoutput,figname)

                 


            
            

            


