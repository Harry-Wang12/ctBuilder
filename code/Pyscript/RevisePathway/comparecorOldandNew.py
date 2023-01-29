import os
from scipy import stats
import statistics


def compareoriginalandnew(alldir):
    inputdir=alldir+"RevisedScoringanaly/"
    summayfilepath = alldir+"OldandNewsummary_withreference.txt"
    randomscorefile= alldir+"randomScoresAnaly/"
    with open(summayfilepath,"w") as summaryresult:
        g = os.walk(inputdir)
        for path,dir_list,file_list in g:  
            for file_name in file_list:
                print(file_name)
                referencelist = {}
                with open(randomscorefile+file_name,"r") as referencefile:
                    for line in referencefile.readlines():
                        if len(line)>1:
                            linedata = line.strip().split("\t")
                            print(linedata[1])
                            p2name = linedata[0]
                            p1name = linedata[1]
                            del(linedata[0])
                            del(linedata[0])
                            del(linedata[0])
                            if p2name not in referencelist.keys():
                                referencelist[p2name]={}
                            referencelist[p2name][p1name]={
                                "mean":0,
                                "value":[]
                            }
                            referencelist[p2name][p1name]["mean"]=0

                            for stringvalue in linedata:
                                referencelist[p2name][p1name]["value"].append(float(stringvalue))
                            referencelist[p2name][p1name]["mean"]=statistics.mean(referencelist[p2name][p1name]["value"])

                with open(inputdir+file_name,"r") as analysisfile:
                   
                    summayanalyonefile={
                        "newaddlargest":{},
                        "newaddAll":[],
                        "original":[],
                        "p2name":"",
                        "filename":file_name,
                        "details":""                      
                    }
                    isnewaddedhighest=0

                    for line in analysisfile.readlines():
                        linedata = line.strip().split("\t")
                        p2name = linedata[0]
                        p1name = linedata[1]
                        summayanalyonefile["details"]= linedata[2].split("~")[0]
                        p1sourceID = p1name.split(",")[0].split(":")[1]
                        if "_x" in p1sourceID:
                            isnewadded=True
                        else:
                            isnewadded=False                        

                        allgenesrelated = linedata[2].split("~")[5]
                        summayanalyonefile["p2name"] = p2name
                        if linedata[3]=="nan":
                            coorelation ="nan"
                            pvalue = "nan"
                        else:
                            coorelation  =float(linedata[3])
                            ttestresult = stats.ttest_1samp(referencelist[p2name][p1name]["value"],coorelation)
                            pvalue = ttestresult.pvalue
                        
                        

                        # index = 0
                        # if coorelation< referencelist[p2name][p1name]["mean"]:
                        #     for value in referencelist[p2name][p1name]["value"]:
                        #         if value<coorelation:
                        #             index+=1
                        # if coorelation > referencelist[p2name][p1name]["mean"]:
                        #     for value in referencelist[p2name][p1name]["value"]:
                        #         if value>coorelation:
                        #             index+=1

                        # pvalue = index/len(referencelist[p2name][p1name]["value"])




                        # kde = stats.gaussian_kde(referencelist[p2name][p1name], bw_method=0.001)

                        # pvalue = kde.integrate_box_1d(coorelation, 1.1)



                        # if linedata[4]=="nan":
                        #     pvalue =1
                        # else:
                        #     pvalue = linedata[4]

                        
                        if isnewadded:
                            addlist={}
                            addlist["p1name"] = p1name
                            addlist["coorelation"] = coorelation
                            addlist["pvalue"] = pvalue
                            addlist["allgene"] = allgenesrelated
                            summayanalyonefile["newaddAll"].append(addlist)
                            if abs(coorelation)>=isnewaddedhighest:
                                summayanalyonefile["newaddlargest"]["p1name"] = p1name
                                summayanalyonefile["newaddlargest"]["coorelation"] = coorelation
                                summayanalyonefile["newaddlargest"]["pvalue"] = pvalue
                                summayanalyonefile["newaddlargest"]["allgene"] = allgenesrelated
                                isnewaddedhighest = abs(coorelation)

                        else:
                            originallist ={}
                            originallist["p1name"] = p1name
                            originallist["coorelation"] = coorelation
                            originallist["pvalue"] = pvalue
                            originallist["allgene"] = allgenesrelated
                            summayanalyonefile["original"].append(originallist)

                    # 
                    detaillist = analysisdetail(summayanalyonefile["details"])
                    pathwaynameoffic = summayanalyonefile["filename"].split("_")[1].split("]")[1]

                    

                    writeline = summayanalyonefile["filename"]+"\n"+summayanalyonefile["details"]+"\n"+pathwaynameoffic+"\n"+summayanalyonefile["p2name"]+"\n"
                    if "p1name" in summayanalyonefile["newaddlargest"].keys():
                        xindex = summayanalyonefile["newaddlargest"]["p1name"].split(",")[0].split("_")[1]
                        writeline += summayanalyonefile["newaddlargest"]["p1name"]+"\t"+summayanalyonefile["newaddlargest"]["allgene"] +"\t"+ str(summayanalyonefile["newaddlargest"]["coorelation"])+"\t"+str(summayanalyonefile["newaddlargest"]["pvalue"])+"\t"+detaillist[xindex]+"\n"
                    if len(summayanalyonefile["original"]) >0:
                        for original in summayanalyonefile["original"]:
                            writeline += original["p1name"]+"\t"+original["allgene"] +"\t"+ str(original["coorelation"])+"\t"+str(original["pvalue"])+"\t"+pathwaynameoffic+"\n"
                    if len(summayanalyonefile["newaddAll"]) >0:
                        for original in summayanalyonefile["newaddAll"]:
                            xindex = original["p1name"].split(",")[0].split("_")[1]
                            writeline += original["p1name"]+"\t"+original["allgene"] +"\t"+ str(original["coorelation"])+"\t"+str(original["pvalue"])+"\t"+detaillist[xindex]+"\n"
                    summaryresult.write(writeline+'\n')


def analysisdetail(detailstring):
    returnlist={}
    
    stringlist=detailstring.strip().split(":")
    for i in range(1,len(stringlist)):
        returnlist["x"+str(i)] = stringlist[i].split("_")[0].split(".")[0]
    
    return returnlist




                        
if __name__=="__main__":

    alldir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_inflamtory/"

    compareoriginalandnew(alldir)

    # alldir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-13-2021_NASH_rescale/"

    # compareoriginalandnew(alldir)
    
    # inputdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-11-2021_Normal_rescale/RevisedScoringanaly/"
    # compareoriginalandnew(inputdir)




