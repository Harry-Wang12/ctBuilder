
import os

def stickrandomresult(folderlist,outputdir):

    if not os.path.exists(outputdir):
        os.makedirs(outputdir) 


    g = os.walk(folderlist[0])

    for path,dir_list,file_list in g:  
        for file_name in file_list:
            print(file_name)  
            featurelist1=[]
            featurelist2=[]
            samplelists=[]
            for loadfolder in folderlist:
                readfilename = loadfolder+file_name
                with open(readfilename,"r") as readfile:
                    filegenelist1=[]
                    filegenelist2=[]
                    isheader = True
                    filesamleplist={}
                    headerorder=[]
                    for line in readfile.readlines():
                        if isheader:
                            linedata = line.strip().split("\t")
                            del(linedata[0])
                            del(linedata[0])
                            for cellname in linedata:
                                filesamleplist[cellname]={}
                                headerorder.append(cellname)
                            isheader=False
                        else:
                            linedata = line.strip().split("\t")
                            genename1 = linedata[0]
                            genename2 = linedata[1]
                            filegenelist1.append(genename1)
                            filegenelist2.append(genename2)
                            del(linedata[0])
                            del(linedata[0])
                            for i in range(len(linedata)):
                                cellname = headerorder[i]
                                filesamleplist[cellname][genename1]=linedata[i]

                    if len(featurelist1)==0:
                        featurelist1=filegenelist1
                        featurelist2=filegenelist2
                    else:
                        featurelist1 = [i for i in featurelist1 if i in filegenelist1]
                        featurelist2 = [i for i in featurelist2 if i in filegenelist2] 

                    for cellname,featurevalue in filesamleplist.items():
                        samplelists.append(featurevalue)

            outputfile = outputdir+file_name
            with open(outputfile,"w") as outputfile:
                title = "RouteID\tRoute"
                for i in range(len(samplelists)):
                    title+="\tCell"+str(i)
                
                outputfile.write(title+"\n")

                for i in range(len(featurelist1)):
                    pline = featurelist1[i]+"\t"+featurelist2[i]

                    for sample in samplelists:
                        pline+="\t"+sample[filegenelist1[i]]
                    
                    outputfile.write(pline+"\n")



            

if __name__=="__main__":



    folderlist=[
        "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/revisedScoring_counts_0_log2/",
        "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/revisedScoring_counts_1_log2/",
        "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/revisedScoring_counts_2_log2/",
        "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/revisedScoring_counts_3_log2/",
        "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/revisedScoring_counts_4_log2/"
    ]
    outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/randomScores/"
    stickrandomresult(folderlist,outputdir)