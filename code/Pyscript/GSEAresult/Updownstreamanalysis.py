import itertools
import os

def loadgraphfromfile(filepath,TFcol=1,target=2,relationcolu=3):
    
    graphdict = {}
    with open(filepath,"r") as graphfile:
        for line in graphfile.readlines():
            linedata= line.strip().split("\t")
            tfgenename = linedata[TFcol-1]
            targetgenename = linedata[target-1]
            relationship = linedata[relationcolu-1]
            if not tfgenename in graphdict.keys():
                graphdict[tfgenename]={}
                    
            if not targetgenename in graphdict[tfgenename].keys():
                graphdict[tfgenename][targetgenename]=[relationship]
            
            elif not relationship in graphdict[tfgenename][targetgenename]:
                graphdict[tfgenename][targetgenename].append(relationship)

    refinegraphdict = {}

    
    for tfgene,targetgene in graphdict.items():
        
        if not tfgene in refinegraphdict.keys():
                refinegraphdict[tfgene]={}
        for targetname,relationships in targetgene.items():

            # if tfgene=='ABL1' and targetname=='CSF1':
            #     print(graphdict['ABL1']['CSF1'])
            # if unkown
            if len(relationships)==3:
                refinegraphdict[tfgene][targetname]=["Activation","Repression"]
            elif len(relationships)==2:
                if "Unknown" in relationships:
                    relationships.remove("Unknown")
                    refinegraphdict[tfgene][targetname]=relationships
                else:
                    refinegraphdict[tfgene][targetname]=["Activation","Repression"]
            elif len(relationships)==1:
                if "Unknown" in relationships:
                    refinegraphdict[tfgene][targetname]=["Activation","Repression"]
                else:
                    refinegraphdict[tfgene][targetname]=relationships

    
            
    return  refinegraphdict



def loadratiofile(filepath):
    CPlist={}
    with open(filepath,'r') as inputfile:
        istitle = True
        for line in inputfile.readlines():
            if not line.startswith("!"):
                if istitle:
                    istitle=False
                else:
                    linedata=line.strip().split('\t')  
                    if len(linedata)>1:                  
                        genename=linedata[0]
                        ratio  = float(linedata[1])
                        CPlist[genename]={
                            'ratio':ratio,
                        }
    return CPlist



def scorelink(tf,target,relation):

    if tf ==0 or target ==0:
        return "bad"


    if relation=="Activation":
        if (tf>0 and target>0) or (tf<0 and target<0):
            return 1
        
    elif relation=="Repression":
        if (tf<0 and target>0) or (tf>0 and target<0):
            return 1

    return 0


# def scoredmap(ratiodict, graphdict):
#     scoredmap = {}

#     for tfgene,targetgene in graphdict.items():
#         if not tfgene in scoredmap.keys():
#                 scoredmap[tfgene]={}
#         for targetname,relationships in targetgene.items():
#             scoredmap[tfgene][targetname]={}
#             # print(tfgene)
#             # print(targetname)
#             # print(relationships)
#             for relationship in relationships:
#                 if tfgene in ratiodict.keys():
#                     tf = ratiodict[tfgene]['ratio']
#                 else:
#                     tf=0

#                 if targetname in ratiodict.keys():
#                     target=ratiodict[targetname]['ratio']
#                 else:
#                     target=0    
                
#                 scoredmap[tfgene][targetname][relationship]=scorelink(tf,target,relationship)

#     return scoredmap






def scoreroute(routepath,graphdict,CPratio,upperstreamscore):

    if float(upperstreamscore)>0:
        signal = 1
    else:
        signal = -1

    finalroutes=[]   


    # get expectation
    # nowsign = signal

    for i in range(len(routepath)-1):

        tf = routepath[i]
        target = routepath[i+1]
        addedlist = []
        if "Activation" in  graphdict[tf][target]:
            addedlist.append({
                'tf':tf,
                'target': target,
                'relation':"Activation",
                
            })


        if "Repression" in  graphdict[tf][target]:
            addedlist.append({
                'tf':tf,
                'target': target,
                'relation':"Repression",
                
            })

        finalroutes.append(addedlist)


    combinations =  list(itertools.product(*finalroutes))
    returnpath = {}

    # start calculate

    for combination in combinations:
        scores =[]
        pathroute= []
        name=""

        nowsign = signal

        for i in range(len(combination)):
            
            if i ==0:
                pathroute.append(combination[i]['tf'])
                pathroute.append(combination[i]['relation'])
                pathroute.append(combination[i]['target'])


            else:
                pathroute.append(combination[i]['relation'])
                pathroute.append(combination[i]['target'])
            name+=combination[i]['relation']



            if combination[i]['relation']=="Activation":
                needsign = nowsign
            elif combination[i]['relation']=="Repression":
                needsign = -1*nowsign

            if combination[i]['target'] in CPratio.keys():
                if CPratio[combination[i]['target']]['ratio'] * needsign>0:
                    scores.append(1)
                else:
                    scores.append(0)
            
            needsign = nowsign

        if len(scores)>0:        
            returnpath[name] = {
                'name':name,
                'pathwaydetails':pathroute,
                'score':(sum(scores)/len(scores) * signal + float(upperstreamscore))/2
            }
        else:
            returnpath[name] = {
                'name':name,
                'pathwaydetails':pathroute,
                'score':0
            }     



    

    return returnpath
        
def loadroutes(routefile):
    paths = {}
    with open(routefile,'r') as routeinput:
        for line in routeinput.readlines():
            linedata = line.strip().split("\t")
            pathname = linedata[0]
            del(linedata[0])
            del(linedata[0])
            if len(linedata)>2:
                paths[pathname] = linedata
    
    return paths


def getPSscore(PSfilepath):
    PSlist = {}
    filenamelist = []
    with open(PSfilepath,'r') as PSfile:
        istitle = True
        for line in PSfile.readlines():
            line = line.replace("\"","")
            if istitle:
                istitle = False
                linedata= line.strip().split('\t')
                del(linedata[0])
                del(linedata[0])
                filenamelist=linedata
            else:
                linedata=line.strip().split('\t')
                kgmlname =linedata[0]
                
                if not kgmlname in PSlist.keys():
                    PSlist[kgmlname] = []
                PSitems = {
                    "genelist" : linedata[1].strip().split(","),
                    "PSscore": {}
                }

                del(linedata[0])
                del(linedata[0])

                for i in range(len(linedata)):
                    PSitems["PSscore"][filenamelist[i]] = linedata[i]
                
                PSlist[kgmlname].append(PSitems)

    return PSlist




def writepathtofile(pathnamedetailslist, outputfilepath,scoredict,filenamelist):
    with open(outputfilepath,'w') as outputfile:
        title = "pathwayname\tpathwaydetails"
        for filename in filenamelist:
            title+='\t'+filename
        outputfile.write(title+'\tSUM\n')

        for pathwayname,pathwaydetails in pathnamedetailslist.items():
            pline = pathwayname+'\t'+', '.join(pathwaydetails)
            sum = 0
            for score in scoredict[pathwayname]:
                pline+='\t'+str(score)
                sum+=score

            outputfile.write(pline+'\t'+str(sum)+'\n')






if __name__=="__main__":

    
    

    graphpath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/CoexpressionDB/TRRUST/trrust_rawdata.human.tsv"

    graphdict = loadgraphfromfile(graphpath,TFcol=1,target=2,relationcolu=3)

    ratiofiepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/testHFD_renew/selected_combineindivpaper/"
    # ratiofiepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/Random_data/random1000_log2ratio_indiv/"


    routepathdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/betweenroute_summary/"
    # routepath = "C:/Users/how17003/Documents/codes/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_7-24/uniquepaths.txt"

    # routes = loadroutes(routepath)

    outputdir= "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/Scoredroutes_withtreatpaper/"
    # outputdir= "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-3/Scoredroutes_random/"


    PSfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/PSresult_cohortresult_small_withtreatpaper.txt"
    # PSfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-3/PSresult_random.txt"
    # routepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_7-29/betweenroute_NASH_RELATED_10_allpath.txt"
    # outputfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_7-24/uniquepathconsistent_random.txt"
    # outputfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_7-29/uniquepathconsistent_sign.txt"

    PSlist = getPSscore(PSfilepath)

    for pathwaypath, pathway_dirs,pathway_files in os.walk(routepathdir):
        
        for pahtway_file in pathway_files:
            print(pahtway_file)
            summaryscore ={}
            routedetails = {}
            filenamelist=[]
            routepath = routepathdir+pahtway_file
            outputfile = outputdir+pahtway_file
            routes = loadroutes(routepath)    
            for path,dir_list,file_list in os.walk(ratiofiepath) :  
                for filename in file_list:  
                    # print(filename)
                    filenamelist.append(filename)
                    CPratio = loadratiofile(ratiofiepath+filename)
                    # scoredmapdict = scoredmap(CPratio, graphdict) 
                    for routename, routepath in routes.items():

                        PSroutescores = PSlist[pahtway_file.strip().split("_")[0]]
                        # PSroutescores = PSlist[pahtway_file.strip()]
                        index = 0
                        for pathwayscores in  PSroutescores:                            
                            if routepath[0] in pathwayscores["genelist"]: 
                                upperstreamscore =  pathwayscores["PSscore"][filename]      
                                # upperstreamscore =  pathwayscores["PSscore"][filename.split(".")[0]]          
                                scoredroutes = scoreroute(routepath,graphdict,CPratio,upperstreamscore)

                                for pathsubname,pathdetails in scoredroutes.items():   
                                    if routename+pathsubname+str(index) not in summaryscore.keys():
                                        summaryscore[routename+pathsubname+str(index)]=[pathdetails['score']]
                                        routedetails[routename+pathsubname+str(index)]=pathdetails['pathwaydetails']
                                    else:
                                        summaryscore[routename+pathsubname+str(index)].append(pathdetails['score'])

                            
                            index +=1
            writepathtofile(routedetails, outputfile,summaryscore,filenamelist)