import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()


def analysisresult(filepath,outputfile,pathwaylistfile,outputpathwaylist,outputfileindex,focuslistfile,threshold):
    heatmatrix = []
    allpathwaylist=[]

    with open(pathwaylistfile,"r") as pathwaylist:
        for line in pathwaylist.readlines():
            pathwayname = line.strip().split(" - ")[1].split(".")[0]
            allpathwaylist.append(pathwayname)

    with open(filepath,"r") as matrixfile:
        index = 0
        for line in matrixfile.readlines():
            line = line.strip() 
            linedata =line.strip().split("\t")
            if len(line)>0:                
                index+=1
                if index==3:
                    P2pathway = line.strip()
                    if P2pathway not in allpathwaylist:  
                        allpathwaylist.append(P2pathway)
                if len(linedata)==5:                   
                    P1pathway = linedata[4]
                    combinepathway = [P2pathway,P1pathway]
                    correlation = abs(float(linedata[2]))
                    if P1pathway not in allpathwaylist:  
                        allpathwaylist.append(P1pathway)
                    
                    needadd = True
                    for combineitem in heatmatrix:
                        if sorted(combineitem["pathwayconbine"]) == sorted(combinepathway):
                            needadd=False
                            if correlation>=combineitem["cor"]:
                                combineitem["cor"]=correlation                                
                                break
                    if needadd:
                        addinitem={
                            "pathwayconbine":sorted([P2pathway,P1pathway]),
                            "cor":correlation
                        }
                        heatmatrix.append(addinitem)


                    # if P1pathway not in heatmatrix[P2pathway].keys():
                    #     heatmatrix[P2pathway][P1pathway]=correlation
                    # else:
                    #     if correlation>=heatmatrix[P2pathway][P1pathway]:
                    #         heatmatrix[P2pathway][P1pathway]=correlation
                    
            if len(line)==0:
                index=0

    allpathwaylist = sorted(allpathwaylist)

    
    with open(focuslistfile,"w") as focuspathway:
        focuslist=[]
        with open(outputfile,"w") as outputmatrix:       
                pline ="Pathway"
                for pathway in allpathwaylist:
                    pline+="\t"+pathway
                pline+="\n"
                for pathway1 in allpathwaylist:
                    pline+=pathway1
                    for pathway2 in allpathwaylist:
                        writefocuse = True
                        checkelist = [pathway1,pathway2]
                        coor = 0
                        for combineitem in heatmatrix:
                            if sorted(combineitem["pathwayconbine"]) == sorted(checkelist):
                                coor=combineitem["cor"]
                        if coor >= threshold:
                            for combineitem in focuslist:
                                if sorted(combineitem["pathwayconbine"]) == sorted(checkelist):
                                    writefocuse = False
                        if coor >= threshold and writefocuse:
                            focuslist.append( {
                                "pathwayconbine":sorted(checkelist),
                                "cor":coor
                            })
                        pline+="\t"+str(coor)
                    pline+="\n"
                outputmatrix.write(pline)

        for combineitem in focuslist:
            focuspathway.write(combineitem["pathwayconbine"][0]+"\t"+ combineitem["pathwayconbine"][1]+"\t"+str(combineitem["cor"])+"\n")


    

    with open(outputfileindex,"w") as outputmatrixindex:       
            pline ="Pathway"
            for i in  range(len(allpathwaylist)):
                pline+="\t"+str(i+1)
            pline+="\n"
            for i in  range(len(allpathwaylist)):
                pline+=str(i+1)
                for j in  range(len(allpathwaylist)):
                    checkelist = [allpathwaylist[i],allpathwaylist[j]]
                    coor = 0
                    for combineitem in heatmatrix:
                        if sorted(combineitem["pathwayconbine"]) == sorted(checkelist):
                            coor=combineitem["cor"]
                    pline+="\t"+str(coor)


                pline+="\n"
            outputmatrixindex.write(pline)

    with open(outputpathwaylist,"w") as outpathwaylist:
        for i in range(len(allpathwaylist)):
            outpathwaylist.write(str(i)+"\t"+allpathwaylist[i]+"\n")
            




    
def drawheatmap(heatmapfile,savefigpath):
    df= pd.read_csv(heatmapfile , delimiter = "\t", header=0,index_col=0)
    fig, ax = plt.subplots(figsize=(20,10))  
    sns.heatmap(df,square=True,xticklabels=False,)
    # plt.imshow(df,cmap='hot',interpolation='nearest')
    plt.savefig(savefigpath)
    # plt.show()



if __name__ =="__main__":
    pathwaylistfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/pathways_list.txt"
    threshold=0.75
    


    # filepathdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/"
    filepathdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_inflamtory/"
    outputpathwaylist =  filepathdir+"heatmappathwaylist.txt"
    filepath = filepathdir+"OldandNewsummary_withreference.txt"
    outputfile = filepathdir+"heatmapmatrix.txt"
    outputfileindex = filepathdir+"heatmapmatrixindex.txt"
    savefigpath = filepathdir+"heatmap.jpg"
    focuslistfile = filepathdir+"focuslist.txt"
    

    analysisresult(filepath,outputfile,pathwaylistfile,outputpathwaylist,outputfileindex,focuslistfile,threshold)

    drawheatmap(outputfileindex,savefigpath)




            













