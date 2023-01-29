

def valueminus(filepath1,filepath2,outputfilepath):
    filelist1={}
    filelist2={}
    with open(filepath1,'r') as file1:
        for line in file1.readlines():
            linedata= line.strip().split('\t')
            target = linedata[0]
            TF = linedata[1]
            importance = float(linedata[2])
            if not target in filelist1.keys():
                filelist1[target]={
                    TF:importance
                }
            else:
                filelist1[target][TF]=importance
    

    with open(filepath2,'r') as file2:
        for line in file2.readlines():
            linedata= line.strip().split('\t')
            target = linedata[0]
            TF = linedata[1]
            importance = float(linedata[2])
            if not target in filelist2.keys():
                filelist2[target]={
                    TF:importance
                }
            else:
                filelist2[target][TF]=importance
    
    with open(outputfilepath,'w')as outputfile:
        
        for target,TFs in filelist1.items():
            for TF,importance in TFs.items():
                if target in filelist2.keys() and TF in filelist2[target].keys():
                    outputfile.write(target+'\t'+TF+'\t'+str(filelist1[target][TF]-filelist2[target][TF])+'\n')
    

def finedelink(filelists,outputfilepath,threshold=0.5):
    
    realdiffdict={}
    diffdict={}


    for filepath in filelists:
        with open(filepath,'r') as inputfile:
            for line in inputfile.readlines():
                linedata = line.strip().split('\t')
                target = linedata[0]
                tf = linedata[1]
                importancediff = float(linedata[2])
                if not target in diffdict.keys():
                    diffdict[target]={}
                if not tf in  diffdict[target].keys():
                    diffdict[target][tf]=[importancediff]
                else:
                    diffdict[target][tf].append(importancediff)

    with open(outputfilepath,'w') as outputfile:
        for target,tfs in diffdict.items():
            for tf,importancelists in tfs.items():
                # importancelists is a list
                plusimp=0
                minusimp = 0
                for value in importancelists:
                    if value >=threshold:
                        plusimp+=1
                    if value <= -1*threshold:
                        minusimp+=1
                outputfile.write(target+'\t'+tf+'\t'+str(plusimp/len(importancelists))+'\t'+str(minusimp/len(importancelists))+"\t"+str(max([plusimp/len(importancelists),minusimp/len(importancelists)]))+'\n')
               












if  __name__=="__main__":

    outputdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_6-27/"

    cresultweightfile1 = outputdir+"kmeans_GSE169445_control_DC1.txt"
    tresultweightfile1 = outputdir+"kmeans_GSE169445_NASH_DC1.txt"
    
    cresultweightfile2 = outputdir+"kmeans_GSE169445_control_DC2.txt"
    tresultweightfile2 = outputdir+"kmeans_GSE169445_NASH_DC2.txt"
    
    cresultweightfile3 = outputdir+"kmeans_GSE169445_control_DC3.txt"
    tresultweightfile3 = outputdir+"kmeans_GSE169445_NASH_DC3.txt"

    cresultweightfile4 = outputdir+"kmeans_GSE169445_control_DC4.txt"
    tresultweightfile4 = outputdir+"kmeans_GSE169445_NASH_DC3.txt"
    
    cresultweightfile5 = outputdir+"kmeans_GSE169445_control_DC5.txt"
    tresultweightfile5 = outputdir+"kmeans_GSE169445_NASH_DC5.txt"
    
   





    outputfilepath1=outputdir+"GSE169445_control_NASH_minus_DC1.txt"
    outputfilepath2=outputdir+"GSE169445_control_NASH_minus_DC2.txt"
    outputfilepath3=outputdir+"GSE169445_control_NASH_minus_DC3.txt"
    outputfilepath4=outputdir+"GSE169445_control_NASH_minus_DC4.txt"
    outputfilepath5=outputdir+"GSE169445_control_NASH_minus_DC5.txt"




    # valueminus(cresultweightfile1,tresultweightfile1,outputfilepath1)
    # valueminus(cresultweightfile2,tresultweightfile2,outputfilepath2)
    # valueminus(cresultweightfile3,tresultweightfile1,outputfilepath3)

    # valueminus(cresultweightfile4,tresultweightfile4,outputfilepath4)
    # valueminus(cresultweightfile5,tresultweightfile5,outputfilepath5)


    filelists = [
        outputdir+"GSE169445_control_NASH_minus_DC1.txt",
        outputdir+"GSE169445_control_NASH_minus_DC2.txt",
        outputdir+"GSE169445_control_NASH_minus_DC3.txt",
        outputdir+"GSE169445_control_NASH_minus_DC4.txt",
        outputdir+"GSE169445_control_NASH_minus_DC5.txt"

    ]

    outputfilepath = outputdir+"kmeans_all.txt"
    finedelink(filelists,outputfilepath,threshold=0.5)