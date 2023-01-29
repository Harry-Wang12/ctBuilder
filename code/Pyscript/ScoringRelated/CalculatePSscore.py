

def calculatePS(inputfilepath,outputfilepath,thresholdP):

    with open(outputfilepath,"w") as outputfile:
        with open(inputfilepath,"r") as inputfile:
            istitle = True
            for line in inputfile.readlines():
                if istitle:
                    istitle=False
                else:
                    significat = 0
                    linedata=line.strip().split("\t")
                    routeid = linedata[0]
                    route = linedata[1]
                    del(linedata[0])
                    del(linedata[0])
                    pline=routeid+"\t"+route
                    for value in linedata:
                        if float(value)<= thresholdP:
                            significat+=1
                    pline+="\t"+str(significat/len(linedata))+"\n"
                    outputfile.write(pline)

if __name__=="__main__":


    thresholdP=0.05
    inputfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-11-2021_NASH_rescale/originalroute_pvalues.txt"
    outputfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-11-2021_NASH_rescale/originalroute_PS.txt"
    calculatePS(inputfilepath,outputfilepath,thresholdP)



