


import os
import scipy.stats as st
import matplotlib.pyplot as plt
import statistics  as stats
import numpy as np
import random



def comparetheresult(sourcedir,resultdir,imageprocessdir):
    if not os.path.exists(resultdir):
        os.makedirs(resultdir)  

    if not os.path.exists(imageprocessdir):
        os.makedirs(imageprocessdir)    

    g = os.walk(sourcedir)  

    for path,dir_list,file_list in g:  
        for file_name in file_list:

            # analysis filename:
            filenamedata = file_name.strip().split("_")
            p2source = filenamedata[3]
            p2target = filenamedata[4]

            selectedp2=[]
            selectedp1list=[]
            print(file_name)

            with open(os.path.join(path, file_name),"r") as scoresfile:

                imgdir = os.path.splitext(imageprocessdir+file_name)[0]+"/"
                if not os.path.exists(imgdir):
                    os.makedirs(imgdir)

                title = True

                for line in scoresfile.readlines():
                    if title:
                        title=False
                    else:
                        linedata = line.strip().split("\t")
                        informationdata = linedata[0].strip().split("~")
                        routepart = informationdata[1]
                        source = informationdata[2]
                        target = informationdata[3]
                        if routepart=="p2" and source==p2source and target==p2target:
                            for i in range(2,len(linedata)):
                                selectedp2.append(float(linedata[i]))
                                # create fig

                        if routepart =="p1" and target==p2source:
                            p1items = {
                                'source':source,
                                'target':target,
                                'description':linedata[1],
                                'valuelist':[],
                            }
                            for i in range(2,len(linedata)):
                                p1items['valuelist'].append(float(linedata[i]))

                            selectedp1list.append(p1items)
            
            if  len(selectedp2)>0 and len(selectedp1list)>0:
                with open(os.path.join(resultdir, file_name),"w") as resultfile:
                    p2title = "p2Source: "+p2source+",p2Target: "+p2target
                    printline=""
                    for p1route in selectedp1list:
                        if not stats.stdev(p1route["valuelist"]) ==0 and  not stats.stdev(selectedp2) ==0:
                            p1title = "p1Source: "+p1route["source"]+",p1Target: "+p1route["target"]+"\t"+p1route['description']                            
                            r,p =st.pearsonr(selectedp2, p1route["valuelist"]) 
                            printline+=p2title+"\t"+p1title+"\t"+str(r)+"\t"+str(p)+"\n"
                            # draw image
                            figname = imgdir+"p1Source_"+p1route["source"]+"_p1Target_"+p1route["target"]+".jpg".replace(",","_")
                            # kwargs = dict(histtype='stepfilled', alpha=0.3, bins=50)

                            m, b = np.polyfit(p1route["valuelist"], selectedp2, 1)
                            # m = slope, b=intercept
                            
                            plt.figure()
                            plt.scatter(p1route["valuelist"], selectedp2)
                            plt.plot([-1,1], m*np.array([-1,1]) + b)


                            # plt.hist(selectedp2, **kwargs)
                            # plt.hist(p1route["valuelist"], **kwargs)
                            plt.savefig(figname)
                            plt.close()

                    resultfile.write(printline)



def comparetheresultwithsubset(sourcedir,resultdir,subsetnumber,iteration):
    if not os.path.exists(resultdir):
        os.makedirs(resultdir)  



    g = os.walk(sourcedir)  

    for path,dir_list,file_list in g:  
        for file_name in file_list:

            # analysis filename:
            filenamedata = file_name.strip().split("_")
            p2source = filenamedata[3]
            p2target = filenamedata[4]

            selectedp2=[]
            selectedp1list=[]
            print(file_name)

            with open(os.path.join(path, file_name),"r") as scoresfile:
                
                # if not os.path.exists(imgdir):
                #     os.makedirs(imgdir)

                title = True

                for line in scoresfile.readlines():
                    if title:
                        title=False
                    else:
                        linedata = line.strip().split("\t")
                        informationdata = linedata[0].strip().split("~")
                        routepart = informationdata[1]
                        source = informationdata[2]
                        target = informationdata[3]
                        if routepart=="p2" and source==p2source and target==p2target:
                            for i in range(2,len(linedata)):
                                selectedp2.append(float(linedata[i]))
                                # create fig

                        if routepart =="p1" and target==p2source:
                            p1items = {
                                'source':source,
                                'target':target,
                                'description':linedata[1],
                                'valuelist':[],
                            }
                            for i in range(2,len(linedata)):
                                p1items['valuelist'].append(float(linedata[i]))

                            selectedp1list.append(p1items)
            
            if  len(selectedp2)>0 and len(selectedp1list)>0:
                with open(os.path.join(resultdir, file_name),"w") as resultfile:
                    p2title = "p2Source: "+p2source+",p2Target: "+p2target
                    printline=""
                    for p1route in selectedp1list:
                        if not stats.stdev(p1route["valuelist"]) ==0 and not stats.stdev(selectedp2) ==0:
                            p1title = "p1Source: "+p1route["source"]+",p1Target: "+p1route["target"]+"\t"+p1route['description'] 
                            printline+=p2title+"\t"+p1title 
                            for i in range(iteration):
                                
                                selectedp2sub = []
                                selectedp1sub = []

                                selectedindex =random.sample(range(len(selectedp2)),subsetnumber)

                                for index in selectedindex:
                                    selectedp2sub.append(selectedp2[index])
                                    selectedp1sub.append(p1route["valuelist"][index])
                                if not stats.stdev(selectedp2sub) ==0 and not stats.stdev(selectedp1sub) ==0:
                                    r,p =st.pearsonr(selectedp2sub, selectedp1sub) 
                                    printline+="\t"+str(r)
                                else:
                                    printline+="\t0"
                            printline+="\n"
                        else:
                            p1title = "p1Source: "+p1route["source"]+",p1Target: "+p1route["target"]+"\t"+p1route['description'] 
                            printline+=p2title+"\t"+p1title 
                            for i in range(iteration):
                                printline+="\t0"
                            printline+="\n"
                    resultfile.write(printline+"\n")




                            # # draw image
                            # figname = imgdir+"p1Source_"+p1route["source"]+"_p1Target_"+p1route["target"]+".jpg".replace(",","_")
                            # # kwargs = dict(histtype='stepfilled', alpha=0.3, bins=50)

                            # m, b = np.polyfit(p1route["valuelist"], selectedp2, 1)
                            # # m = slope, b=intercept
                            
                            # plt.figure()
                            # plt.scatter(p1route["valuelist"], selectedp2)
                            # plt.plot([-1,1], m*np.array([-1,1]) + b)


                            # # plt.hist(selectedp2, **kwargs)
                            # # plt.hist(p1route["valuelist"], **kwargs)
                            # plt.savefig(figname)
                            # plt.close()

                    




if __name__=="__main__":

    sourcedir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/randomScores/"
    resultdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_non_inflamtory/randomScoresAnaly/"
    
    # comparetheresult(sourcedir,resultdir,imageprocessdir)
    subsetnumber=50
    iteration=1000

    comparetheresultwithsubset(sourcedir,resultdir,subsetnumber,iteration)

    # sourcedir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-13-2021_NASH_rescale/revisedScoring/"
    # resultdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-13-2021_NASH_rescale/RevisedScoringanaly/"
    # imageprocessdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-13-2021_NASH_rescale/RevisedScoringanaly_image/"
    # comparetheresult(sourcedir,resultdir,imageprocessdir)

    # sourcedir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-11-2021_NASH_rescale/revisedScoring/"
    # resultdir="C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-11-2021_NASH_rescale/RevisedScoringanaly/"
    # comparetheresult(sourcedir,resultdir)






    # print(os.path.join(path, dir_name))
    

    