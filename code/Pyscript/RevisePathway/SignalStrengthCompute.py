
import os
import scipy.stats as st

def getSignalScore(sourcedir,resultdir):
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

                    p2line = p2title

                    for value in selectedp2:
                        p2line+="\t"+str(value)
                    
                    p2line+="\n"
                    resultfile.write(p2line)

                    
                    for p1route in selectedp1list:
                        
                        p1title = "p1Source: "+p1route["source"]+",p1Target: "+p1route["target"]+"\t"+p1route['description']                            
                        r,p =st.pearsonr(selectedp2, p1route["valuelist"]) 
                        printline=p1title

                        for value in p1route["valuelist"]:
                            printline+="\t"+str(value)
                    
                        printline+="\n"
                        resultfile.write(printline)
                                              

                    





    